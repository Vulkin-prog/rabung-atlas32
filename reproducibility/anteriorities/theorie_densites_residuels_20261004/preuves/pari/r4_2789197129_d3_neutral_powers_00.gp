default(parisizemax,512000000);
default(parisize,128000000);
run()={my(fp=(5955787396)*y^0+(1165516800)*y^1+(57532228)*y^2+(-153600)*y^3+(-15167)*y^4+(0)*y^5+(1)*y^6,nf=nfinit(fp),alpha=Mod((-2932530727481/198543389353)*y^0+(52940280149/1985433893530)*y^1+(43539840/198543389353)*y^2+(-57547387/794173557412)*y^3+(-17280/198543389353)*y^4+(22767/3970867787060)*y^5,fp),vs=[[0, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0], [0, 1, 2, 3, 0, 0, 0, 2, 0, 0, 1, 3, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 3, 0, 2, 0, 0, 0, 0, 1]],z,rt,t);
 for(i=1,#vs,z=prod(j=1,#vs[i],(alpha+j-1)^vs[i][j]);
  if(!nfeltispower(nf,z,4,&rt),error("not fourth power"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,fp)^4!=z,error("identity"));print("POWER ",vector(6,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
