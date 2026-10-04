default(parisizemax,512000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(1)*x^4+(126)*x^3+(5593)*x^2+(102996)*x^1+(657004)*x^0,nf=nfinit(P),sets=[[2, 14, 24, 25]],z,rt,t);
 for(i=1,#sets,z=prod(j=1,#sets[i],Mod(x+sets[i][j],P));
  if(!nfeltissquare(nf,z,&rt),error("unsatisfied relation"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,P)^2!=z,error("identity"));print("SQUARE ",vector(4,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
