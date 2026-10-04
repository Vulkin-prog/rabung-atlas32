default(parisizemax,1024000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(4)*x^4+(168)*x^3+(2485)*x^2+(15006)*x^1+(30001)*x^0,Q,C,fp,nf,fac,roots,raw);
 Q=subst(P,x,x/4)*4^(poldegree(P)-1);C=polcompositum(Q,Q,1);
 for(i=1,#C,if(poldegree(C[i][1])==4*(4-1),raw=C[i][1];break));fp=polredbest(raw);print("STAGE pair");
 
 C=polcompositum(fp,x^2+1,1);raw=C[1][1];for(i=2,#C,if(poldegree(C[i][1])>poldegree(raw),raw=C[i][1]));fp=polredbest(raw);print("STAGE cyclotomic ",poldegree(fp));
 fp=subst(fp,x,y);nf=nfinit(fp);print("STAGE nfinit");fac=nffactor(nf,P);roots=[];
 for(i=1,matsize(fac)[1],if(poldegree(fac[i,1])==1,roots=concat(roots,[-polcoef(fac[i,1],0)/polcoef(fac[i,1],1)])));
 if(#roots!=2,error("root number"));
 print("FIELD ",vector(poldegree(fp)+1,i,Str(polcoef(fp,i-1))));
 print("ROOTS ",vector(2,i,vector(poldegree(fp),j,Str(polcoef(lift(roots[i]),j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
