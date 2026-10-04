"""Elementary finite-module algorithms; no PARI dependency."""

def kernel4(rows):
    A=[[x%4 for x in row] for row in rows];s=len(A);m=len(A[0]);T=[[int(i==j) for j in range(m)] for i in range(m)]
    k=0;odd=0
    def swap_columns(a,b):
        for row in A:row[a],row[b]=row[b],row[a]
        for row in T:row[a],row[b]=row[b],row[a]
    for phase in (1,2):
        while k<m and k<s:
            loc=next(((i,j) for i in range(k,s) for j in range(k,m) if (A[i][j]%2==1 if phase==1 else A[i][j]==2)),None)
            if loc is None:break
            i,j=loc;A[i],A[k]=A[k],A[i];swap_columns(j,k)
            if phase==1:
                u=pow(A[k][k],-1,4);A[k]=[(u*x)%4 for x in A[k]]
            for i in range(s):
                if i==k:continue
                c=A[i][k] if phase==1 else int(A[i][k]==2)
                if c:A[i]=[(a-c*b)%4 for a,b in zip(A[i],A[k])]
            for j in range(m):
                if j==k:continue
                c=A[k][j] if phase==1 else int(A[k][j]==2)
                if c:
                    for row in A:row[j]=(row[j]-c*row[k])%4
                    for row in T:row[j]=(row[j]-c*row[k])%4
            k+=1
        if phase==1:odd=k
    even=k-odd;gens=[]
    for j in range(odd,m):
        scale=2 if j<k else 1
        v=[(scale*T[i][j])%4 for i in range(m)]
        assert any(v) and all(sum(a*b for a,b in zip(row,v))%4==0 for row in rows)
        gens.append(v)
    assert all(A[i][j]==(1 if i==j and i<odd else 2 if i==j and i<k else 0) for i in range(s) for j in range(m))
    return {'C4':odd,'C2':even,'log2_order':2*odd+even,'order':4**odd*2**even,'kernel_generators':gens}

def echelon(rows):
    pivots={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in pivots:x^=pivots[p]
            else:pivots[p]=x;break
    return pivots

def nullspace(matrix):
    n=len(matrix[0]);rows=[sum((x%2)<<i for i,x in enumerate(row)) for row in matrix]
    pivots=echelon(rows);free=[i for i in range(n) if i not in pivots];basis=[]
    for i in free:
        v=1<<i
        for p in sorted(pivots):
            if (pivots[p]&v).bit_count()%2:v|=1<<p
        assert all((x&v).bit_count()%2==0 for x in rows)
        basis.append(v)
    return basis

def intersection(A,B):
    A=list(echelon(A).values());B=list(echelon(B).values());pivots={};answers=[]
    for j,original in enumerate(A+B):
        v=original;tag=1<<j
        while v:
            p=v.bit_length()-1
            if p in pivots:w,t=pivots[p];v^=w;tag^=t
            else:pivots[p]=(v,tag);break
        if not v:
            z=0
            for i,a in enumerate(A):
                if tag>>i&1:z^=a
            if z:answers.append(z)
    return list(echelon(answers).values())
