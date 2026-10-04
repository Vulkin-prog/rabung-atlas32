default(parisizemax,512000000);
default(parisize,128000000);
run()={my(fp=(792100)*y^0+(0)*y^1+(94161)*y^2+(0)*y^3+(606)*y^4+(0)*y^5+(1)*y^6,nf=nfinit(fp),alpha=Mod((-2999/324)*y^0+(0)*y^1+(353/1296)*y^2+(0)*y^3+(1/1296)*y^4+(0)*y^5,fp),vs=[[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 3, 1]],z,rt,t);
 for(i=1,#vs,z=prod(j=1,#vs[i],((alpha+j)/alpha)^vs[i][j]);
  if(!nfeltispower(nf,z,4,&rt),error("not fourth power"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,fp)^4!=z,error("identity"));print("POWER ",vector(6,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
