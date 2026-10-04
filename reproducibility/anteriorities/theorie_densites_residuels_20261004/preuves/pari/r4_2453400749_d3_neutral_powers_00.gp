default(parisizemax,512000000);
default(parisize,128000000);
run()={my(fp=(143650)*y^0+(42970)*y^1+(3673)*y^2+(-648)*y^3+(-104)*y^4+(-2)*y^5+(1)*y^6,nf=nfinit(fp),alpha=Mod((-4723/420)*y^0+(-3793/6825)*y^1+(-2897/18200)*y^2+(-211/18200)*y^3+(11/54600)*y^4+(1/7800)*y^5,fp),vs=[[3, 0, 0, 0, 0, 3, 0, 2, 0, 0, 3, 3, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], [1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 3, 0, 2, 1, 0, 3, 0, 0, 0, 0, 0, 3, 0, 3, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 3, 2, 0, 3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [3, 3, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 3, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 1]],z,rt,t);
 for(i=1,#vs,z=prod(j=1,#vs[i],(alpha+j-1)^vs[i][j]);
  if(!nfeltispower(nf,z,4,&rt),error("not fourth power"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,fp)^4!=z,error("identity"));print("POWER ",vector(6,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
