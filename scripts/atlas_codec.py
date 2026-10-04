"""Read-only original atlas codec, extracted without changing its arithmetic."""
import hashlib, struct, zlib
MAGIC=b"VDWSM1\0\0"
HEADER=struct.Struct("<8sQQQQII32s")

def expected_mask(p):return sum(1<<(r-2) for r in range(2,10) if (p-1)%r==0)

def minimum_k(p,r,m,l):return max(m+1,2*l+2 if ((p-1)//2)%r==0 else l+2)

def get_uint(raw,index):
 value=0
 for n in range(10):
  if index>=len(raw):raise ValueError('truncated integer')
  b=raw[index];index+=1;value|=(b&127)<<(7*n)
  if b<128:
   if value>=2**64 or (n and b==0):raise ValueError('noncanonical integer')
   return value,index
 raise ValueError('oversized integer')

def decode_frame(frame):
 if len(frame)<HEADER.size:raise ValueError('truncated header')
 magic,n,size,lo,hi,count,flags,h=HEADER.unpack(frame[:HEADER.size])
 if magic!=MAGIC or flags not in [0,1] or n!=len(frame)-HEADER.size:raise ValueError('frame header')
 if size>128*1024**2 or n>128*1024**2:raise ValueError('oversized frame')
 decomp=zlib.decompressobj();raw=decomp.decompress(frame[HEADER.size:],size+1)
 if len(raw)!=size or not decomp.eof or decomp.unused_data or hashlib.sha256(raw).digest()!=h:raise ValueError('frame hash/length')
 rows=[];i=0;p=lo-1
 for unused in range(count):
  delta,i=get_uint(raw,i);p+=delta
  if delta==0 or not lo<=p<hi or i>=len(raw):raise ValueError('record ordering')
  mask=raw[i];i+=1
  if mask!=expected_mask(p):raise ValueError('decoded mask')
  prof=[]
  for r in range(2,10):
   if mask&(1<<(r-2)):
    m,i=get_uint(raw,i);l,i=get_uint(raw,i);v={'r':r,'max':m,'prefix':l,'K':minimum_k(p,r,m,l)}
    if not 1<=l<=m<=(p-1)//r:raise ValueError('decoded profile bound')
    if flags:
     v['witness'],i=get_uint(raw,i)
     if not 1<=v['witness']<=(p-1)//2:raise ValueError('decoded witness')
    prof.append(v)
  rows.append({'p':p,'profiles':prof})
 if i!=len(raw):raise ValueError('trailing record bytes')
 return rows,{'lo':lo,'hi':hi,'count':count,'witnesses':bool(flags),'raw_sha256':h.hex(),'raw_bytes':size,'compressed_bytes':n}

def science_hash(rows):
 h=hashlib.sha256(b'vdw-p-profile-v1-u64le')
 for row in rows:
  p=row['p'];h.update(struct.pack('<QB',p,expected_mask(p)))
  for v in row['profiles']:h.update(struct.pack('<QQQ',v['max'],v['prefix'],v['K']))
 return h.hexdigest()

FRAME=struct.Struct("<8sQQIIII32s32s32s")
def decode(frame):
    if len(frame)<FRAME.size:raise ValueError('A4 short frame')
    magic,lo,hi,n,csize,rsize,flags,ph,sh,digest=FRAME.unpack(frame[:FRAME.size])
    if magic!=b'VDWA4V1\0' or n>4096 or rsize!=4*n or flags or csize!=len(frame)-FRAME.size:raise ValueError('A4 header')
    d=zlib.decompressobj();raw=d.decompress(frame[FRAME.size:],rsize+1)
    if len(raw)!=rsize or not d.eof or d.unused_data or hashlib.sha256(raw).digest()!=digest:raise ValueError('A4 frame hash')
    return list(struct.unpack('<'+'I'*n,raw)),dict(lo=lo,hi=hi,count=n,prime_hash=ph.hex(),profile_science_hash=sh.hex(),values_sha256=digest.hex())
