default(parisizemax,512000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(1)*x^3+(-197)*x^2+(-9264)*x^1+(-90252)*x^0,nf=nfinit(P),sets=[[2, 12, 15, 17, 18, 20, 28, 29, 37, 39, 44], [3, 12, 18, 22, 23, 24, 26, 27, 37, 43, 46]],z,rt,t);
 for(i=1,#sets,z=prod(j=1,#sets[i],Mod(x+sets[i][j],P));
  if(!nfeltissquare(nf,z,&rt),error("unsatisfied relation"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,P)^2!=z,error("identity"));print("SQUARE ",vector(3,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
