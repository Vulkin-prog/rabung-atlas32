default(parisize,64000000);
run()={my(P=(1)*x^4+(100)*x^3+(3258)*x^2+(34844)*x^1+(26309)*x^0,q,fac,types,seen=List(),factors);
 forprime(q=2,101,if(polcoef(P,poldegree(P))%q==0||poldisc(P)%q==0,next);
  fac=factor(Mod(P,q));types=vecsort(vector(matsize(fac)[1],i,poldegree(fac[i,1])));
  if(setsearch(Set(Vec(seen)),types),next);listput(seen,types);
  factors=vector(matsize(fac)[1],i,vector(poldegree(fac[i,1])+1,j,lift(polcoef(fac[i,1],j-1))));
  print("FACTOR ",[q,types,factors]));print("DISCRIMINANT ",poldisc(P));print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
