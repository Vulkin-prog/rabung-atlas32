default(parisizemax,512000000);
default(parisize,128000000);
run()={my(P=(1)*x^4+(41)*x^3+(424)*x^2+(164)*x^1+(-8320)*x^0,q=1000000,rr,iota,logs,k=0);
 while(k<384,q=nextprime(q+1);if(q%4!=1||(1>1&&q%1==0),next);rr=Vec(vecsort(lift(polrootsmod(P,q))));if(#rr<4,next);rr=rr[1..4];
  if(prod(i=1,4,prod(j=0,28-1,Mod(rr[i]+j,q)))==0,next);
  iota=lift(znprimroot(q)^((q-1)/4));logs=vector(4*28,j,my(t=lift(Mod(rr[(j-1)\28+1]+(j-1)%28,q)^((q-1)/4)));if(t==1,0,if(t==iota,1,if(t==q-1,2,3))));
  print("SAMPLE ",[q,rr,iota,logs]);k++);
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
