"""Read-only CPU extraction/classification of already retained ordinary runs.
No worker, no prime search, no character scan. Fail closed; never alter repo.
"""
import argparse, csv, hashlib, importlib.util, json, os, pathlib, subprocess, sys, time
PIN = {'scripts/ternary_io.py':'dfdad3beeea7880e4471dc30742e1ac286cd578b',
'data/ternary-million/manifest.json':'12a98fdce8277663ed71a04da1be73fb67fb4c0b',
'data/ternary-million/threshold_counts.csv':'d628fe9a87f4ee098b8e76da59de061575ef839d'}
def gitsha(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def require(ok,msg):
 if not ok: raise ValueError(msg)
def main():
 a=argparse.ArgumentParser();a.add_argument('--repo',required=True);a.add_argument('--output-dir',required=True);a.add_argument('--seconds',type=int,default=450);ns=a.parse_args()
 start=time.monotonic();root=pathlib.Path(ns.repo).resolve();out=pathlib.Path(ns.output_dir).resolve()
 require(not out.is_relative_to(root),'output must be outside checkout');require(not out.exists(),'output directory already exists');out.mkdir(parents=True)
 for path,want in PIN.items():require(gitsha((root/path).read_bytes())==want,'input identity '+path)
 spec=importlib.util.spec_from_file_location('pinned_ternary_io',root/'scripts/ternary_io.py');mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 ds=mod.Dataset(root/'data/ternary-million');require(len(ds.blocks)==248,'block count')
 ts=out/'ordinary_runs.tsv';seen=[0]*4;totals=[[0]*10 for _ in range(4)];canonical=[[0]*10 for _ in range(4)];zeros=[[0]*10 for _ in range(4)];fac2=[[0]*10 for _ in range(4)]
 listsig=hashlib.sha256();samples=[]
 with ts.open('x') as f:
  for e in sorted(ds.blocks.values(),key=lambda d:(d['s'],d['start'])):
   for p in ds.read_block(e):
    s=p.s;require(p.index==seen[s],'index coverage');seen[s]+=1
    for j,k in enumerate(range(18,28)):
     n=p.ordinary_count(k);totals[s][j]+=n;canonical[s][j]+=p.canonical_count(k);zeros[s][j]+=n==0;fac2[s][j]+=n*(n-1)
    for x,l,c in p.runs:
     if p.kind((x,l,c))!='ordinary':continue
     line=f'{s}\t{p.p}\t{x}\t{l}\t{c}\n';f.write(line);listsig.update(line.encode())
     if len(samples)<1000 or (p.index%1024==0):samples.append([s,p.p,x,l,c])
   require(time.monotonic()-start<ns.seconds,'CPU deadline')
 require(seen==[250000]*4,'incomplete population')
 with (root/'data/ternary-million/threshold_counts.csv').open() as f:
  for row in csv.DictReader(f):
   if row['stratum']=='all':continue
   s=int(row['stratum']);j=int(row['threshold'])-18
   require(int(row['profiles'])==seen[s],'CSV profiles')
   for key,arr in [('ordinary_runs',totals),('canonical_runs',canonical),('zero_ordinary_profiles',zeros),('ordinary_factorial_moment_sum',fac2)]: require(int(row[key])==arr[s][j],'CSV '+key)
 src=pathlib.Path(__file__).resolve().parent/'classify.cpp';exe=out/'classify'
 subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra',str(src),'-o',str(exe)],check=True,timeout=60)
 subprocess.run([str(exe),str(ts),str(out/'HITS.csv'),str(out/'SUMMARY.csv')],check=True,timeout=max(1,ns.seconds-int(time.monotonic()-start)))
 with (out/'SUMMARY.csv').open() as f:
  for row in csv.DictReader(f):require(int(row['ordinary_runs'])==totals[int(row['stratum'])][int(row['threshold'])-18],'classifier total')
 # Verify every positive witness with direct integer arithmetic, separately from C++.
 hits=0
 with (out/'HITS.csv').open() as f:
  for row in csv.DictReader(f):
   hits+=1;p=int(row['p']);x=int(row['start']);l=int(row['length']);mask=int(row['mask'])
   for prefix,bound,dbound,bit in [('small',64,256,1),('wide',128,4096,2)]:
    if mask&bit:
     j,A,D=(int(row[prefix+'_'+s]) for s in ['j','A','D']);require(0<=j<l and abs(A)<=bound and abs(D)<=dbound,'witness bounds')
     require(D!=0 and (D<0 or __import__('math').isqrt(D)**2!=D),'rational discriminant');require((A*A-D)%4==0,'noninteger polynomial')
     y=x+j;require((y*y+A*y+(A*A-D)//4)%p==0,'quadratic witness')
   if mask&4:
    q=int(row['q']);aa=int(row['a']);require(2<=q<=128 and __import__('math').gcd(aa,q)==1 and q*x<=aa*p<=q*(x+l-1),'rational witness')
 (out/'SAMPLE_RUNS.json').write_text(json.dumps(samples,separators=(',',':'))+'\n')
 # Only compact text artifacts leave the runner. The temporary enumeration is removed.
 size=ts.stat().st_size;ts.unlink();exe.unlink()
 report={'schema':'ternary-coordinate-classification-v1','profiles':sum(seen),'blocks':248,'ordinary_runs21':sum(t[3] for t in totals),'tagged_runs21':hits,'source_git_blobs':PIN,'data_commit':'452907920dc623d4d1e6ae7675428fe02ed840c4','workflow_commit':os.environ.get('GITHUB_SHA'),'run_id':os.environ.get('GITHUB_RUN_ID'),'native_reader_full_validation':True,'does_not_recompute_characters_or_maxima':True,'input_tsv_sha256':listsig.hexdigest(),'input_tsv_bytes':size,'seconds':time.monotonic()-start,'catalogues':{'quadratic_small':{'Amax':64,'Dmax':256},'quadratic_wide':{'Amax':128,'Dmax':4096},'rational_qmax':128},'status':'PASS'}
 (out/'EXPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
 for name in ['classify.cpp','export_classified.py']:(out/name).write_bytes((pathlib.Path(__file__).resolve().parent/name).read_bytes())
 require(time.monotonic()-start<ns.seconds,'CPU deadline');require(sum(p.stat().st_size for p in out.iterdir())<64*1024**2,'output budget')
 (out/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted(out.iterdir())))
if __name__=='__main__':main()
