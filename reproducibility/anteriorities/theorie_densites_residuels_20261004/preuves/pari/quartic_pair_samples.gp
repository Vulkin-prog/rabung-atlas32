default(parisizemax,512000000);
default(parisize,64000000);
run()={my(P=(1)*x^4+(100)*x^3+(3258)*x^2+(34844)*x^1+(26309)*x^0,q=1000000,roots,logs,k=0,v);
 while(k<384,q=nextprime(q+1);roots=Vec(vecsort(lift(polrootsmod(P,q))));if(#roots<2,next);
  roots=roots[1..2];if(prod(i=1,2,prod(j=0,52-1,Mod(roots[i]+j,q)))==0,next);
  logs=vector(2*52,k,my(i=(k-1)\52+1,j=(k-1)%52,v=lift(Mod(roots[i]+j,q)^((q-1)/2)));if(v==1,0,1));
  print("SAMPLE ",[q,roots,logs]);k++);
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
