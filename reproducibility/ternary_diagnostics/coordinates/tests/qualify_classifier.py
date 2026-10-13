import csv, json, math, pathlib, random, subprocess, tempfile
BASE=pathlib.Path(__file__).resolve().parents[1]
def brute(p,x,l,bound,dbound):
 for b in range(-bound,bound+2*(l-1)+1):
  v=pow(2*x+b,2,p)
  D=v if v<=dbound else v-p
  if abs(D)>dbound or D==0 or (D>0 and math.isqrt(D)**2==D) or (b*b-D)%4:continue
  j=max(0,(b-bound+1)//2)
  return (j,b-2*j,D)
 return None
def ref(row):
 s,p,x,l,c=row
 small=brute(p,x,l,64,256);wide=brute(p,x,l,128,4096)
 r=None
 for q in range(2,129):
  a=(q*x+p-1)//p
  if math.gcd(a,q)==1 and q*x<=a*p<=q*(x+l-1):r=(q,a);break
 return (small is not None)+2*(wide is not None)+4*(r is not None),small,wide,r
rng=random.Random(61501);data=[]
for p in [10007,100003,1000003,1000000007,8589934609,17179869043]:
 for i in range(100):
  l=rng.randrange(21,90);x=rng.randrange(1,(p-1)//2-l+1);data.append((0,p,x,l,i%3))
 for q in range(3,100):
  aa=rng.randrange(1,(q+1)//2);x=(aa*p+q-1)//q-5
  if x>=1 and x+25<(p-1)//2:data.append((0,p,x,26,0))
# Engineered roots close to catalogue boundaries.
p=1000000007
for A in [-128,-64,-3,0,5,64,128]:
 for D in [-4095,-255,-31,73,97,4093]:
  if (A*A-D)%4:continue
  sq=pow(D%p,(p+1)//4,p)
  if sq*sq%p!=D%p:continue
  for rt in [sq,p-sq]:
   for j in [0,20,25]:
    x=((rt-A)*pow(2,-1,p)-j)%p
    if 1<=x and x+25<=(p-1)//2:data.append((0,p,x,26,0))
data=sorted(data,key=lambda d:(d[0],d[1],d[2]));exps=[ref(d) for d in data]
with tempfile.TemporaryDirectory() as tmp:
 tmp=pathlib.Path(tmp);exe=tmp/'classify';subprocess.run(['g++','-O2','-std=c++17','-Wall','-Wextra','-fsanitize=undefined',str(BASE/'src/classify.cpp'),'-o',str(exe)],check=True)
 inp=tmp/'input';inp.write_text(''.join('\t'.join(map(str,d))+'\n' for d in data))
 subprocess.run([str(exe),str(inp),str(tmp/'hits'),str(tmp/'summary')],check=True,capture_output=True)
 # Synthetic records can duplicate keys, compare full ordered positive rows.
 rows=list(csv.DictReader((tmp/'hits').open()));expected=[(d,e) for d,e in zip(data,exps) if e[0]]
 assert len(rows)==len(expected)
 for row,(d,e) in zip(rows,expected):
  assert tuple(int(row[k]) for k in ['stratum','p','start','length','colour'])==d
  assert int(row['mask'])==e[0]
  for prefix,v in [('small',e[1]),('wide',e[2])]:
   if v:assert tuple(int(row[prefix+'_'+k]) for k in ['j','A','D'])==v
  if e[3]:assert (int(row['q']),int(row['a']))==e[3]
 # Include a very long synthetic interval to exercise the former 32-bit square risk.
 inp.write_text('0\t1000000007\t12345678\t65535\t0\n');subprocess.run([str(exe),str(inp),str(tmp/'hits'),str(tmp/'summary')],check=True,capture_output=True)
 for bad in ['0 8191 1 21 0\n','0 10007 5000 21 0\n','4 10007 1 21 0\n','0 17179869184 1 21 0\n','garbage\n']:
  inp.write_text(bad);assert subprocess.run([str(exe),str(inp),str(tmp/'hits'),str(tmp/'summary')],capture_output=True).returncode!=0
report={'random_and_engineered_intervals':len(data),'positive_intervals':len(expected),'max_length_65535_ubsan':'PASS','negative_tests':5,'status':'PASS'}
(BASE/'results/CLASSIFIER_TESTS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
