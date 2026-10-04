default(parisize,128000000);default(parisizemax,512000000);default(realprecision,64);setrand(1);
run()={my(P=x^4+80*x^3+2281*x^2+27240*x+116640,fp,bn,roots,fac,z,t,qs=[],S,U,E,units,h,xx,
 q=1000000,rr,logs,k=0,v,nsamples=224);
 gettime();fp=subst(polredbest(polcompositum(P,x^2-265,1)[1][1]),x,y);
 bn=bnfinit(fp,1);fac=nffactor(bn,P);z=-polcoef(fac[1,1],0);
 roots=[z,-40-z];for(i=1,4,t=-polcoef(fac[i,1],0);if(t!=roots[1]&&t!=roots[2],roots=concat(roots,[t,-40-t]);break));
 print("FIELD ",vector(poldegree(fp)+1,j,Str(polcoef(fp,j-1))));
 print("ROOTS ",vector(#roots,i,vector(poldegree(fp),j,Str(polcoef(lift(roots[i]),j-1)))));
 for(j=-11,51,qs=concat(qs,Vec(factor(abs(subst(P,x,-j)))[,1])));qs=Set(qs);
 S=concat(vector(#qs,j,idealprimedec(bn,qs[j])));U=bnfunits(bn,S);h=#U[1];
 E=vector(#roots,i,Mat(vector(52,j,bnfisunit(bn,roots[i]+j-1,U))));
 units=vector(h,i,nffactorback(bn,U[1][i]));
 for(i=1,#roots,if(matsize(E[i])[1]!=h,error("decomposition shape"));for(j=1,52,
  if(nffactorback(bn,units,E[i][,j])!=nfalgtobasis(bn,roots[i]+j-1),error("explicit identity"))));
 if(nfeltpow(bn,units[h],bn.tu[1])!=nfalgtobasis(bn,1),error("torsion"));
 print("TORSION ",bn.tu[1]);
 print("MATRICES ",vector(#roots,k,vector(h,i,vector(52,j,E[k][i,j]))));
 print("UNITS ",vector(h,i,xx=lift(nfbasistoalg(bn,units[i]));vector(poldegree(fp),j,Str(polcoef(xx,j-1)))));
 print("EXPLICIT_IDENTITIES 1");
 while(k<nsamples,q=nextprime(q+1);if(q>=2^32,error("range"));
  rr=lift(polrootsmod(fp,q));if(!#rr,next);xx=rr[1];
  if(prod(i=1,#roots,prod(j=0,51,Mod(subst(lift(roots[i]),y,xx)+j,q)))==0,next);
  logs=vector(#roots*52,j,my(i=(j-1)\52+1,tj=(j-1)%52,v=lift(Mod(subst(lift(roots[i]),y,xx)+tj,q)^((q-1)/2)));if(v==1,0,1));
  print("SAMPLE ",[q,xx,logs]);k++);
 print("MILLISECONDS ",gettime());print("END 1");
};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
