default(parisizemax,512000000);
default(parisize,64000000);
run()={my(P=(1)*x^3+(-197)*x^2+(-9264)*x^1+(-90252)*x^0,q=1000000,roots,logs,k=0,v);
 while(k<384,q=nextprime(q+1);roots=Vec(vecsort(lift(polrootsmod(P,q))));if(#roots<3,next);
  roots=roots[1..3];if(prod(i=1,3,prod(j=0,53-1,Mod(roots[i]+j,q)))==0,next);
  logs=vector(3*53,k,my(i=(k-1)\53+1,j=(k-1)%53,v=lift(Mod(roots[i]+j,q)^((q-1)/2)));if(v==1,0,1));
  print("SAMPLE ",[q,roots,logs]);k++);
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
