import json, math, sys

def factor(n):
    out=[];d=2
    while d*d<=n:
        if n%d==0:
            out.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:out.append(n)
    return out

def certify(row):
    p=row['p'];s=row['twist_parameter'];assert 1<=s<p
    assert s==1 or pow(s,(p-1)//2,p)==p-1
    def add(a,b):
        if a is None:return b
        if b is None:return a
        x,y=a;u,v=b
        if x==u:
            if (y+v)%p==0:return None
            assert y==v;slope=(3*x*x-20*s*x+9*s*s)*pow(2*y,-1,p)%p
        else:slope=(v-y)*pow(u-x,-1,p)%p
        z=(slope*slope+10*s-x-u)%p
        return z,(slope*(x-z)-y)%p
    def multiply(a,n):
        result=None
        while n:
            if n&1:result=add(result,a)
            a=add(a,a);n>>=1
        return result
    lcm=1
    for w in row['certificate_points']:
        x,y,n=w['x'],w['y'],w['order'];assert 0<=x<p and 0<=y<p and n>=1
        assert y*y%p==x*(x-s)*(x-9*s)%p
        assert factor(n)==w['prime_factors'];assert multiply((x,y),n) is None
        assert all(multiply((x,y),n//q) is not None for q in w['prime_factors'])
        lcm=math.lcm(lcm,n)
    assert lcm==row['order_lcm'];R=math.isqrt(4*p);lower=p+1-R;upper=p+1+R
    assert row['hasse_lower']==lower and row['hasse_upper']==upper
    first=((lower+lcm-1)//lcm)*lcm;assert first<=upper and first+lcm>upper
    N=first if s==1 else 2*p+2-first;assert row['elliptic_points']==N
    def chi(x):return 1 if pow(x%p,(p-1)//2,p)==1 else -1
    boundary=4*chi(2)+2*chi(6)+2*chi(3)+2*chi(-1)+2*chi(-2)
    numerator=p-10-boundary+N-p-2;assert numerator>=0 and numerator%8==0
    assert row['quadruples']==numerator//8
    return dict(status='PASS',point_certificates=len(row['certificate_points']),twisted_curve=s!=1,order_lcm=lcm)


if __name__ == '__main__':
    data=json.load(open(sys.argv[1]))
    if isinstance(data,dict) and 'results' in data:
        rows=[r.get('elliptic',r) for r in data['results']]
    else:
        rows=data if isinstance(data,list) else [data]
    for row in rows:
        p=row['p']
        assert 1000<=p<2**33 and p%2 and all(p%d for d in range(3,math.isqrt(p)+1,2))
        certify(row)
    print(json.dumps({'status':'PASS','certificates':len(rows)}))
