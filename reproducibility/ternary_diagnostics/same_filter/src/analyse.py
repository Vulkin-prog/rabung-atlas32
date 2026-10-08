from pathlib import Path
from fractions import Fraction as F
import csv,json,hashlib
from filter_model import parameters,mean_upper,innovation_bounds,coordinate_mask,filtered_mean,is_prime_trial,require
ROOT=Path(__file__).resolve().parents[1]


def enc(x):
    if isinstance(x,F):return {'fraction':str(x),'decimal':float(x)}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def main():
    path=ROOT/'input/COMPACT.json'
    require(hashlib.sha256(path.read_bytes()).hexdigest()=='27f9679d7b2cdbb3884ef5a201092a43c1a1e992b4d2b5900c8dc30915389e48','compact input changed')
    compact=json.loads(path.read_text())
    hits=[{k:int(v) for k,v in row.items()} for row in csv.DictReader((ROOT/'input/HITS.csv').open())]
    require(len(hits)==2300 and len({(x['p'],x['start']) for x in hits})==2300,'hit list shape')
    require(len({x['p'] for x in hits})==2187,'hit prime count')
    par=parameters();small=parameters(64,256,128)
    catalogs={'small':(small['L'],small['D'],1),'quadratic':(par['L'],par['D'],2),
              'rational':(0,par['R'],4),'union':(par['L'],par['J'],6)}
    totals={};lam={};marked={};n=250000
    for s in range(4):
        st=compact['strata'][s];require(st['profiles']==n,'stratum size')
        for k in range(21,28):
            totals[s,k]=sum(int(x)*count for x,count in st['histograms'][str(k)].items())
            lam[s,k]=F(st['sum_p']-n*(2*k+4),3**k)
            marked[s,k]=sum(x['stratum']==s and x['length']>=k and bool(x['mask']&6) for x in hits)
    S={(int(x['stratum']),int(x['threshold']),x['category']):x for x in csv.DictReader((ROOT/'input/SUMMARY.csv').open())}
    for s in range(4):
        for k in range(21,28):require(marked[s,k]==int(S[s,k,'union_wide']['tagged_runs']),'summary mismatch')
    envelopes=[]
    for name,(L,J,bit) in catalogs.items():
        for k in range(21,28):
            observed=sum(x['length']>=k and bool(x['mask']&bit) for x in hits)
            bound=1000000*mean_upper(k,3,L,J)
            envelopes.append({'catalog':name,'k':k,'observed_marked_runs':observed,
                              'iid_expected_sum_upper':bound,
                              'observed_minus_iid_lower':F(observed)-bound,
                              'observed_over_iid_lower':F(observed)/bound})
    contrasts=[];strata=[]
    for k in range(21,28):
        ub=500000*mean_upper(k,3,par['L'],par['J'])
        raw=sum((1 if s%2==0 else -1)*(totals[s,k]-lam[s,k]) for s in range(4))
        removed=sum((1 if s%2==0 else -1)*marked[s,k] for s in range(4))
        contrasts.append({'k':k,'initial_contrast':raw,'tagged_raw_contrast':removed,
                          'background_contrast_absolute_upper':ub,
                          'tagged_centered_contrast':(F(removed)-ub,F(removed)+ub),
                          'residual_centered_contrast':(raw-removed-ub,raw-removed+ub),
                          'fraction_of_initial_centered_contrast':((removed-ub)/raw,(removed+ub)/raw)})
        for s in range(4):
            u=ub/2
            old=totals[s,k]-marked[s,k]-lam[s,k]
            strata.append({'s':s,'k':k,'count':totals[s,k],'marked':marked[s,k],
                           'iid_whole_mean_without_terminal':lam[s,k],
                           'remaining_count':totals[s,k]-marked[s,k],
                           'remaining_excess_same_filter_bounds':(old,old+u)})
    continuation=[]
    for a in range(21,27):
        for b in range(a+1,28):
            q=F(1,3**(b-a));raw=F(0);tag=F(0)
            for s in range(4):
                sign=1 if s%2==0 else -1
                raw+=sign*(totals[s,b]-q*totals[s,a]+F(2*(b-a)*n,3**b))
                tag+=sign*(marked[s,b]-q*marked[s,a])
            lo,hi=innovation_bounds(a,b,3,par['J']);margin=500000*(hi-lo)
            coarse=500000*(mean_upper(b,3,par['L'],par['J'])+q*mean_upper(a,3,par['L'],par['J']))
            row={'a':a,'b':b,'initial_contrast':raw,'tagged_raw_contrast':tag,
                 'iid_tagged_innovation_each_profile_bounds':(lo,hi),
                 'iid_contrast_absolute_upper':margin,'coarse_bound':coarse,
                 'tagged_centered_contrast':(tag-margin,tag+margin),
                 'residual_centered_contrast':(raw-tag-margin,raw-tag+margin)}
            if raw>0:row['centered_fraction']=[(tag-margin)/raw,(tag+margin)/raw]
            continuation.append(row)
    samples=sorted(set([next(x['p'] for x in hits if x['stratum']==s) for s in range(4)]+[8642031919]))
    evaluations=[]
    for p in samples:
        require(is_prime_trial(p),'sample prime trial division failed')
        mask=coordinate_mask(p)
        vals=[]
        for k in [21,24,27]:
            lo,hi=filtered_mean(mask['h'],k,3,mask['components'])
            bound=mean_upper(k,3,par['L'],par['J'])
            require(0<=lo<=hi<=bound,'sample model outside uniform envelope')
            vals.append({'k':k,'iid_marked_mean_bounds':(lo,hi),'width':hi-lo,
                         'uniform_upper':bound})
        evaluations.append({'p':p,'not_representative':True,'split_D':mask['split_discriminants'],
                            'input_intervals':len(mask['intervals']),'rational_cuts':len(mask['cut_edges']),
                            'complement_components':len(mask['components']),
                            'ordinary_domain_length':mask['h']-2,
                            'unmarked_vertices':sum(mask['components']),'values':vals})
    result={'status':'PASS','data_commit':'452907920dc623d4d1e6ae7675428fe02ed840c4',
            'scope':'Same geometric filter, independent uniform three-colour words; original prime strata fixed. No data acquisition; no causal or prospective claim.',
            'terminal_ordinary_mean_error_sum_upper':F(1,10**100),
            'terminal_note':'Rows use (p-2k-4)/3^k; add an outward error <10^-100 for the omitted positive iid terminal over the sample.',
            'parameters':par,'small_parameters':small,'envelopes':envelopes,
            'contrasts':contrasts,'stratum_residuals':strata,'continuations':continuation,
            'sample_exact_geometry_evaluations':evaluations}
    out=ROOT/'results/ANALYSE.json';out.write_text(json.dumps(enc(result),indent=2,sort_keys=True)+'\n')
    print(json.dumps(enc({'parameters':par,'main_envelopes':[x for x in envelopes if x['catalog']=='union' and x['k'] in [21,24,27]],
                          'passages':[x for x in continuation if (x['a'],x['b']) in [(21,24),(24,27)]],
                          'sample_p':samples}),indent=2))

if __name__=='__main__':main()
