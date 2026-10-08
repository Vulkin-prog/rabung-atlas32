import collections,csv,hashlib,json,math,pathlib
from fractions import Fraction as F
from relations import certificate,verify_certificate
BASE=pathlib.Path(__file__).resolve().parents[1]
def enc(x):
    if isinstance(x,F):return {'fraction':str(x),'decimal':float(x)}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x

def main():
    raw=(BASE/'input/COMPACT.json').read_bytes();assert hashlib.sha256(raw).hexdigest()=='27f9679d7b2cdbb3884ef5a201092a43c1a1e992b4d2b5900c8dc30915389e48'
    compact=json.loads(raw);folder=BASE/'input/coordinate_export'
    for line in (folder/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ');assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest
    export=json.loads((folder/'EXPORT.json').read_text());assert export['status']=='PASS' and export['profiles']==1000000 and export['ordinary_runs21']==1234808
    summary=list(csv.DictReader((folder/'SUMMARY.csv').open()));S={(int(r['stratum']),int(r['threshold']),r['category']):{k:int(v) if k!='category' else v for k,v in r.items()} for r in summary}
    hits=[{k:int(v) for k,v in r.items()} for r in csv.DictReader((folder/'HITS.csv').open())];assert len(hits)==2300
    assert len({(r['stratum'],r['p'],r['start']) for r in hits})==len(hits)
    bit={'quadratic_small':1,'quadratic_wide':2,'rational':4,'union_wide':6}
    for s in range(4):
        for k in range(21,28):
            a=sum(int(x)*n for x,n in compact['strata'][s]['histograms'][str(k)].items())
            for cat,b in bit.items():
                rr=[r for r in hits if r['stratum']==s and r['length']>=k and r['mask']&b]
                assert S[s,k,cat]['ordinary_runs']==a and S[s,k,cat]['tagged_runs']==len(rr)
                assert S[s,k,cat]['affected_profiles']==len({r['p'] for r in rr})
    models={};annot=[]
    for row in hits:
        if not row['mask']&2:continue
        j,a,d=row['wide_j'],row['wide_A'],row['wide_D'];u=a+2*j;v=j*j+a*j+(a*a-d)//4;key=f'{u}:{v}:{row["length"]}'
        if key not in models:
            models[key]=certificate(u,v,row['length']);verify_certificate(models[key])
        c=models[key]
        assert c['phase_proof'] is None or row['colour']==0
        assert c['two_proof'] is None or row['stratum']%2==0
        annot.append({'stratum':row['stratum'],'p':row['p'],'start':row['start'],'length':row['length'],'model':key,'neutral_forced':c['phase_proof'] is not None,'two_cube_forced':c['two_proof'] is not None,'neutral_observed':row['colour']==0})
    rows=[];lam={};A={}
    for s in range(4):
        st=compact['strata'][s];n=st['profiles']
        for k in range(18,28):
            A[s,k]=sum(int(x)*n for x,n in st['histograms'][str(k)].items());lam[s,k]=F(st['sum_p']-n*(2*k+4),3**k)
        for k in range(21,28):
            h=S[s,k,'union_wide']['tagged_runs'];rows.append({'stratum':s,'k':k,'count':A[s,k],'reference_sum':lam[s,k],'excess':A[s,k]-lam[s,k],'theta':F(A[s,k])/lam[s,k],'removed':h,'residual_count':A[s,k]-h,'residual_excess':A[s,k]-h-lam[s,k],'residual_theta':F(A[s,k]-h)/lam[s,k],'new_voids':S[s,k,'union_wide']['fully_covered_profiles']})
    contrasts=[]
    for k in range(21,28):
        cube=sum(A[s,k]-lam[s,k] for s in [0,2]);non=sum(A[s,k]-lam[s,k] for s in [1,3]);hc=sum(S[s,k,'union_wide']['tagged_runs'] for s in [0,2]);hn=sum(S[s,k,'union_wide']['tagged_runs'] for s in [1,3]);raw=cube-non
        contrasts.append({'k':k,'cube_excess':cube,'noncube_excess':non,'contrast_count_equivalent':raw,'removed_cube':hc,'removed_noncube':hn,'contrast_after_removal':raw-hc+hn,'fraction_of_contrast_displaced':F(hc-hn)/raw})
    continuations=[]
    for a in range(21,27):
        for b in range(a+1,28):
            q=F(1,3**(b-a));count=removed=F(0)
            by=[]
            for s in range(4):
                raw=F(A[s,b])-q*A[s,a]+F(2*(b-a)*250000,3**b)
                tagged=F(S[s,b,'union_wide']['tagged_runs'])-q*S[s,a,'union_wide']['tagged_runs']
                count+=(-1 if s%2 else 1)*raw;removed+=(-1 if s%2 else 1)*tagged
                by.append({'stratum':s,'innovation_sum':raw,'tagged_sum':tagged,'residual_sum':raw-tagged})
            continuations.append({'a':a,'b':b,'raw_contrast':count,'tagged_contrast':removed,'residual_contrast':count-removed,'by_stratum':by})
    certificates_summary=[]
    for k in range(21,28):
        aa=[r for r in annot if r['length']>=k]
        certificates_summary.append({'k':k,'quadratic_hits':len(aa),'neutral_observed':sum(r['neutral_observed'] for r in aa),'neutral_forced':sum(r['neutral_forced'] for r in aa),'two_cube_forced':sum(r['two_cube_forced'] for r in aa),'both_forced':sum(r['neutral_forced'] and r['two_cube_forced'] for r in aa),'cubic_stratum_hits':sum(r['stratum']%2==0 for r in aa)})
    result={'status':'PASS','scope':'Descriptive finite sample; no causal attribution, no predictive refit. Reference sum omits positive iid terminal term 3^(1-h).','rows':rows,'contrasts':contrasts,'continuations':continuations,'certificates':certificates_summary,'distinct_models':len(models),'distinct_polynomials':len({(c['u'],c['v']) for c in models.values()}),'stored_identities':sum(len(c['relations']) for c in models.values()),'source_export':export}
    (BASE/'results/ANALYSE.json').write_text(json.dumps(enc(result),indent=2,sort_keys=True)+'\n')
    (BASE/'results/CERTIFICATS.json').write_text(json.dumps({'models':models,'occurrences':annot},indent=2,sort_keys=True)+'\n')
    print(json.dumps(enc({'models':len(models),'identities':result['stored_identities'],'certificates':certificates_summary,'contrasts':contrasts,'passages':[v for v in continuations if (v['a'],v['b']) in [(21,24),(24,27)]]}),indent=2))
if __name__=='__main__':main()
