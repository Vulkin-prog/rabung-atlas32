"""Independent local checks of the new export; not a second million-run scan."""
import collections,copy,csv,hashlib,json,math,pathlib,subprocess,tempfile
BASE=pathlib.Path(__file__).resolve().parents[1];D=BASE/'input/coordinate_export'
rows=[{k:int(v) for k,v in row.items()} for row in csv.DictReader((D/'HITS.csv').open())]
summary=list(csv.DictReader((D/'SUMMARY.csv').open()))
for line in (D/'SHA256SUMS').read_text().splitlines():
 digest,name=line.split('  ');assert hashlib.sha256((D/name).read_bytes()).hexdigest()==digest
for name in ['classify.cpp','export_classified.py']:assert (D/name).read_bytes()==(BASE/'src'/name).read_bytes()
# Trial division by an independently built list of all possible prime divisors.
limit=math.isqrt(max(r['p'] for r in rows));sieve=bytearray(b'\1')*(limit+1);sieve[:2]=b'\0\0'
for q in range(2,math.isqrt(limit)+1):
 if sieve[q]:sieve[q*q::q]=b'\0'*(((limit-q*q)//q)+1)
primes=[q for q in range(2,limit+1) if sieve[q]];prime_ids=sorted({r['p'] for r in rows})
for p in prime_ids:
 assert p%6==1
 for q in primes:
  if q*q>p:break
  assert p%q
characters=0;quad=rat=0
for r in rows:
 p,x,l=r['p'],r['start'],r['length'];e=(p-1)//3
 val=pow(x,e,p);characters+=1
 assert (val==1)==(r['colour']==0)
 for j in range(1,l):assert pow(x+j,e,p)==val;characters+=1
 assert pow(x-1,e,p)!=val and pow(x+l,e,p)!=val;characters+=2
 assert (pow(2,e,p)==1)==(r['stratum']%2==0);characters+=1
 for label,ab,db,bit in [('small',64,256,1),('wide',128,4096,2)]:
  if r['mask']&bit:
   j,a,d=(r[label+'_'+z] for z in ['j','A','D']);assert 0<=j<l and abs(a)<=ab and abs(d)<=db
   assert d!=0 and (d<0 or math.isqrt(d)**2!=d) and (a*a-d)%4==0
   assert ((x+j)**2+a*(x+j)+(a*a-d)//4)%p==0;quad+=1
 if r['mask']&4:
  a,q=r['a'],r['q'];assert math.gcd(a,q)==1 and 2<=q<=128 and q*x<=a*p<=q*(x+l-1);rat+=1
# Full direct reclassification of every positive and all exported sample intervals.
def simple_quad(p,x,l,ab,db):
 for t in range(-ab,ab+2*(l-1)+1):
  v=((2*x+t)**2)%p;d=v if v<=db else v-p
  if abs(d)>db or d==0 or (d>0 and math.isqrt(d)**2==d) or (t*t-d)%4:continue
  return True
 return False
def mask_of(s,p,x,l,c):
 bits=int(simple_quad(p,x,l,64,256))+2*int(simple_quad(p,x,l,128,4096))
 for q in range(2,129):
  a=(q*x+p-1)//p
  if q*x<=a*p<=q*(x+l-1) and math.gcd(a,q)==1:bits|=4;break
 return bits
samples=json.loads((D/'SAMPLE_RUNS.json').read_text());seen={(r['stratum'],r['p'],r['start']):r['mask'] for r in rows}
for item in samples:
 assert mask_of(*item)==seen.get(tuple(item[:3]),0)
for r in rows:assert mask_of(*(r[k] for k in ['stratum','p','start','length','colour']))==r['mask']
# Polynomial identities checked by convolution and monic polynomial division,
# independently of the two-coordinate ring multiplication used by generation.
def polyprod(indices):
 f=[1]
 for j in indices:
  g=[0]*(len(f)+1)
  for i,c in enumerate(f):g[i]+=j*c;g[i+1]+=c
  f=g
 return f
def subtract(left,right):
 out=[0]*max(len(left),len(right))
 for i,a in enumerate(left):out[i]+=a
 for i,a in enumerate(right):out[i]-=a
 return out
def identity_ok(rel,u,v):
 f=subtract([rel['den']*a for a in polyprod(rel['left'])],[rel['num']*a for a in polyprod(rel['right'])])
 for i in range(len(f)-1,1,-1):
  a=f[i];f[i]=0;f[i-1]-=a*u;f[i-2]-=a*v
 return all(a==0 for a in f)
def factor(a):
 a=abs(a);out={};q=2
 assert a
 while q*q<=a:
  while a%q==0:out[q]=out.get(q,0)+1;a//=q
  q+=1
 if a>1:out[a]=out.get(a,0)+1
 return out
def relation_vector(rel):
 d={0:len(rel['left'])-len(rel['right'])}
 for val,sign in [(rel['num'],-1),(rel['den'],1)]:
  for p,e in factor(val).items():d[p]=d.get(p,0)+sign*e
 return {p:e%3 for p,e in d.items() if e%3}
def verify_proof(c,name,key):
 proof=c[name]
 if proof is None:return
 out=collections.Counter()
 for i,e in proof:
  assert e in [1,2] and 0<=i<len(c['relations'])
  for p,a in relation_vector(c['relations'][i]).items():out[p]+=e*a
 assert {p:a%3 for p,a in out.items() if a%3}=={key:1}
certs=json.loads((BASE/'results/CERTIFICATS.json').read_text());ids=proofs=0
for c in certs['models'].values():
 assert c['D']==c['u']**2-4*c['v']
 for rel in c['relations']:
  assert rel['num'] and rel['den'] and all(0<=j<c['length'] for j in rel['left']+rel['right'])
  assert identity_ok(rel,c['u'],c['v']);ids+=1
 for name,key in [('phase_proof',0),('two_proof',2)]:
  verify_proof(c,name,key);proofs+=c[name] is not None
# Direct exhaustive A,j search checks the reduction of two loops to one t-loop.
for s,p,x,l,c in samples[:30]:
 for ab,db in [(64,256),(128,4096)]:
  exists=False
  for j in range(l):
   for a in range(-ab,ab+1):
    val=((2*(x+j)+a)**2)%p;d=val if val<=db else val-p
    if abs(d)<=db and d!=0 and (d<0 or math.isqrt(d)**2!=d) and (a*a-d)%4==0:exists=True
  assert exists==simple_quad(p,x,l,ab,db)
# Two explanatory identities printed in the report.
assert identity_ok({'left':[8,18],'right':[6],'num':1,'den':1},25,138)
assert identity_ok({'left':[8,15],'right':[9],'num':-2,'den':1},25,138)
# Counterexamples / deliberate corruption.
c=next(v for v in certs['models'].values() if v['relations'] and v['phase_proof'])
r=copy.deepcopy(c['relations'][0]);r['num']+=1;assert not identity_ok(r,c['u'],c['v'])
bad=copy.deepcopy(c);idx,coef=bad['phase_proof'][0];bad['phase_proof'][0][1]=3-coef
try:verify_proof(bad,'phase_proof',0)
except AssertionError:pass
else:raise AssertionError('mutated proof accepted')
# Recompute integer totals from witness lists, including overlap correction.
for r in summary:
 s,k,cat=int(r['stratum']),int(r['threshold']),r['category'];bit={'quadratic_small':1,'quadratic_wide':2,'rational':4,'union_wide':6}[cat]
 hit=[v for v in rows if v['stratum']==s and v['length']>=k and v['mask']&bit]
 assert len(hit)==int(r['tagged_runs']);assert len({v['p'] for v in hit})==int(r['affected_profiles'])
report={'status':'PASS','positive_runs':len(rows),'prime_trial_divisions_certified':len(prime_ids),'natural_character_evaluations':characters,'quadratic_witness_checks_including_nested_catalogues':quad,'rational_witnesses':rat,'sample_intervals_reclassified':len(samples),'all_positive_intervals_reclassified':len(rows),'nested_parameter_comparisons':60,'polynomial_identities':ids,'linear_implication_certificates':proofs,'distinct_models':len(certs['models']),'summary_rows_reconciled':len(summary),'negative_algebra_tests':2,'does_not_recompute_global_maxima':True,'sample_not_representative':True}
(BASE/'results/VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
