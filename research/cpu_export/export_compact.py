"""Read-only export of existing native ternary data; no scanner is started.

Pinned reader/manifest/CSV identities are mandatory. The complete native reader
validates prime lists, parent frames, counts and positions. Output is a NEW JSON
outside the repository, written only after all 1,000,000 profiles reconcile.
Execution history is recorded separately; no scan or model fitting is performed.
"""
from __future__ import annotations
from collections import Counter
import argparse,csv,hashlib,importlib.util,json,os,signal,sys,tempfile,time
from pathlib import Path
sys.dont_write_bytecode=True
from fractions import Fraction

def require(condition, message):
    if not condition:
        raise ValueError(message)

def cyclic_gaps(n, points):
    points = sorted(points)
    require(n >= 1 and len(points) == len(set(points)), 'invalid circle or repeated point')
    require(all(type(x) is int and 0 <= x < n for x in points), 'point outside circle')
    return tuple([b-a for a,b in zip(points, points[1:])] + [points[0]+n-points[-1]]) if points else ()

def hard_core_Q(n, m, k):
    require(n >= 1 and m >= 0 and k >= 1 and n >= m*k, 'invalid separated configuration')
    if not m:
        return Fraction(0)
    u = n-m*k
    return Fraction(1,m) + Fraction(u*(m-1)*(u+m), m*(m+1)*n*n)

PIN={
 'scripts/ternary_io.py':'dfdad3beeea7880e4471dc30742e1ac286cd578b',
 'data/ternary-million/manifest.json':'12a98fdce8277663ed71a04da1be73fb67fb4c0b',
 'data/ternary-million/threshold_counts.csv':'d628fe9a87f4ee098b8e76da59de061575ef839d'}
K=tuple(range(18,28));GAP_K=(21,24,27);GRID=4096


def git_blob(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()


class Compensated:
    def __init__(self):self.s=0.0;self.c=0.0
    def add(self,value):
        y=value-self.c;t=self.s+y;self.c=(t-self.s)-y;self.s=t


class Reducer:
    def __init__(self):
        self.n=[0]*4
        self.hist=[[Counter() for k in K] for s in range(4)]
        self.joint=[[ [Counter() for b in K] for a in K] for s in range(4)]
        self.sum_p=[0]*4;self.sum_p2=[0]*4;self.sum_p_count=[[0]*10 for s in range(4)]
        self.canonical=[[0]*10 for s in range(4)]
        self.gaps={ (s,k):{'histogram_bins':Counter(),'full_circle_gaps':0,'empty_profiles':0,
                         'pooled_gaps':0,'Q_sum':Compensated(),'Q_reference_sum':Compensated()}
                   for s in range(4) for k in GAP_K}

    def add(self,s,index,p,counts,L,C,ordinary_runs):
        require(type(s)is int and 0<=s<4 and index==self.n[s], 'stratum/order mismatch')
        require(type(p)is int and p>4 and p%2==1, 'odd modulus required')
        counts=tuple(counts);n=(p-1)//2
        require(len(counts)==10 and all(type(v)is int and v>=0 for v in counts)
                and all(a>=b for a,b in zip(counts,counts[1:])), 'invalid count vector')
        runs=sorted(ordinary_runs)
        require(all(type(x)is int and type(Ln)is int and 2<=x and Ln>=21 and x+Ln-1<n for x,Ln in runs),
                'ordinary run not strictly interior or not retained at 21')
        require(all(x+Ln<=y for (x,Ln),(y,Lm) in zip(runs,runs[1:])), 'overlapping runs')
        require(all(sum(Ln>=k for x,Ln in runs)==counts[k-18] for k in range(21,28)),
                'positions disagree with the retained counts')
        self.n[s]+=1;self.sum_p[s]+=p;self.sum_p2[s]+=p*p
        for a,x in enumerate(counts):
            self.hist[s][a][x]+=1
            self.sum_p_count[s][a]+=p*x
            self.canonical[s][a]+=x+int(L>=K[a])+int(C>=K[a])
            for b in range(a+1,len(K)):
                self.joint[s][a][b][(x,counts[b])]+=1
        for k in GAP_K:
            g=cyclic_gaps(n,[x-1 for x,Ln in runs if Ln>=k]);entry=self.gaps[s,k]
            require(not g or min(g)>=k, 'cyclic hard-core separation failed')
            entry['empty_profiles']+=int(not g);entry['pooled_gaps']+=len(g)
            entry['full_circle_gaps']+=int(len(g)==1)
            for d in g:
                if d<n:entry['histogram_bins'][(GRID*d)//n]+=1
            entry['Q_sum'].add(sum((d/n)**2 for d in g))
            entry['Q_reference_sum'].add(float(hard_core_Q(n,len(g),k)))

    def export(self):
        strata=[]
        for s in range(4):
            joints={f'{K[a]}:{K[b]}':{f'{x}:{y}':v for (x,y),v in sorted(self.joint[s][a][b].items())}
                    for a in range(10) for b in range(a+1,10)}
            geometry={}
            for k in GAP_K:
                item=self.gaps[s,k]
                bins=item['histogram_bins']
                lo=sum(v*b*b for b,v in bins.items())+item['full_circle_gaps']*GRID*GRID
                hi=sum(v*(b+1)*(b+1) for b,v in bins.items())+item['full_circle_gaps']*GRID*GRID
                geometry[str(k)]={
                    'non_atom_bins':dict(sorted(bins.items())),
                    **{key:item[key] for key in ('full_circle_gaps','empty_profiles','pooled_gaps')},
                    'Q_sum_numerical':item['Q_sum'].s,
                    'Q_reference_sum_numerical':item['Q_reference_sum'].s,
                    'Q_sum_exact_grid_bounds':{'lower_numerator':lo,'upper_numerator':hi,'denominator':GRID*GRID}}
            strata.append({'stratum':s,'profiles':self.n[s],'sum_p':self.sum_p[s],'sum_p2':self.sum_p2[s],
                           'sum_p_times_count':dict(zip(map(str,K),self.sum_p_count[s])),
                           'histograms':{str(k):dict(sorted(self.hist[s][a].items())) for a,k in enumerate(K)},
                           'canonical_totals':dict(zip(map(str,K),self.canonical[s])),
                           'joint_histograms':joints,'gap_geometry':geometry})
        return {'schema':'ternary-compact-gaps-and-counts-v1','strata':strata,
                'gap_thresholds':GAP_K,'gap_bin_denominator':GRID,
                'gap_weighting':'One gap per ordinary start; empty profiles retained separately; circular closure of half-period.',
                'precision':'Histograms/weighted integer sums/grid envelopes exact; Q sums are compensated floating-point reductions.',
                'source_git_blobs':PIN,'does_not_recompute_maxima_or_characters':True}

    def reconcile(self,csv_data):
        rows=list(csv.DictReader(csv_data.decode().splitlines()))
        for s in range(4):
            require(self.n[s]==250000,'incomplete stratum')
        for row in rows:
            k=int(row['threshold']);a=k-18
            ss=range(4) if row['stratum']=='all' else [int(row['stratum'])]
            hist=Counter()
            for s in ss:hist.update(self.hist[s][a])
            require(sum(hist.values())==int(row['profiles']) and hist[0]==int(row['zero_ordinary_profiles'])
                    and sum(x*v for x,v in hist.items())==int(row['ordinary_runs'])
                    and sum(x*(x-1)*v for x,v in hist.items())==int(row['ordinary_factorial_moment_sum'])
                    and sum(self.canonical[s][a] for s in ss)==int(row['canonical_runs']),
                    'archived CSV reconciliation failed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=int,default=450)
    args=parser.parse_args();repo=args.repo.resolve();output=args.output.absolute()
    require(output.parent.is_dir(),'output parent must already exist')
    require(not output.exists() and not output.is_symlink(),'output already exists')
    require(not output.resolve().is_relative_to(repo),'output must be outside source repository')
    require(1<=args.seconds<=570,'time limit must be in 1..570 seconds')
    def expired(*_):raise TimeoutError('bounded export expired; no completed output accepted')
    signal.signal(signal.SIGALRM,expired);signal.alarm(args.seconds)
    raw={}
    for path,expected in PIN.items():
        raw[path]=(repo/path).read_bytes()
        require(git_blob(raw[path])==expected,'pinned source differs: '+path)
    spec=importlib.util.spec_from_file_location('pinned_ternary_io',repo/'scripts/ternary_io.py')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    data=module.Dataset(repo/'data/ternary-million');reducer=Reducer();start=time.monotonic()
    for profile in data.iter_profiles():
        ordinary=[(x,L) for x,L,c in profile.runs if profile.kind((x,L,c))=='ordinary']
        reducer.add(profile.s,profile.index,profile.p,
                    [profile.ordinary_count(k) for k in K],profile.L,profile.fields[2],ordinary)
        if (profile.index+1)%32768 == 0 or profile.index == 249999:
            print(json.dumps({'stratum':profile.s,'profiles_read':sum(reducer.n),'seconds':round(time.monotonic()-start,3)}), file=sys.stderr, flush=True)
    reducer.reconcile(raw['data/ternary-million/threshold_counts.csv'])
    result=reducer.export();result['elapsed_seconds']=time.monotonic()-start
    result['execution']={'run_id':os.environ.get('GITHUB_RUN_ID'),'head_sha':os.environ.get('GITHUB_SHA'),'data_commit':'452907920dc623d4d1e6ae7675428fe02ed840c4','native_reader_full_validation':True,'blocks_expected':len(data.blocks)}
    payload=(json.dumps(result,sort_keys=True,ensure_ascii=False,separators=(',',':'))+'\n').encode()
    temp=None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent,prefix='.vdw-export-',delete=False) as handle:
            temp=Path(handle.name);handle.write(payload);handle.flush();os.fsync(handle.fileno())
        os.link(temp,output)
    finally:
        if temp is not None:temp.unlink(missing_ok=True)
    signal.alarm(0)
    print(json.dumps({'status':'PASS','profiles':sum(reducer.n),'output':str(output),
                      'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}))
if __name__=='__main__':main()
