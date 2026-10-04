default(parisizemax,512000000);
default(parisize,128000000);
run()={my(P=(5)*x^3+(255)*x^2+(2818)*x^1+(1848)*x^0,q=1000000,rr,iota,logs,k=0);
 while(k<384,q=nextprime(q+1);if(q%4!=1||(5>1&&q%5==0),next);rr=Vec(vecsort(lift(polrootsmod(P,q))));if(#rr<3,next);rr=rr[1..3];
  if(prod(i=1,3,prod(j=0,29-1,Mod(rr[i]+j,q)))==0,next);
  iota=lift(znprimroot(q)^((q-1)/4));logs=vector(3*29,j,my(t=lift(Mod(rr[(j-1)\29+1]+(j-1)%29,q)^((q-1)/4)));if(t==1,0,if(t==iota,1,if(t==q-1,2,3))));
  print("SAMPLE ",[q,rr,iota,logs]);k++);
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
