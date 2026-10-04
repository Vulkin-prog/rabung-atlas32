default(parisize,128000000);setrand(1);
run()={my(P=(1)*x^3+(37)*x^2+(390)*x^1+(1143)*x^0,r=2,m=55,q=1000000,roots,z,t,logs,k=0,bad);
 while(k<192,
 q=nextprime(q+1);if(q>=2^32,error("range"));if((q-1)%r,next);
 roots=lift(polrootsmod(P,q));if(#roots!=3,next);bad=0;
 for(i=1,3,for(j=0,m-1,if((roots[i]+j)%q==0,bad=1)));if(bad,next);
 z=0;for(a=2,q-1,t=lift(Mod(a,q)^((q-1)/r));if(t^r%q==1 && (r!=4 || t^2%q!=1) && t!=1,z=t;break));
 if(!z,error("primitive unity root"));
 logs=vector(3*m,a,my(i=1+(a-1)\m,j=(a-1)%m,v=lift(Mod(roots[i]+j,q)^((q-1)/r)),ell=-1);for(e=0,r-1,if(z^e%q==v,ell=e;break));if(ell<0,error("character logarithm"));ell);
 print("SAMPLE ",[q,Vec(roots),z,logs]);k++;
 );
};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
