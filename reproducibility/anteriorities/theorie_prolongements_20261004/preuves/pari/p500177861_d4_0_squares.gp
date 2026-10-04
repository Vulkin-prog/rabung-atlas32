default(parisizemax,512000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(1)*x^4+(107)*x^3+(4112)*x^2+(66262)*x^1+(379356)*x^0,nf=nfinit(P),sets=[[2, 15, 23, 29, 42, 47], [6, 15, 16, 18, 21, 23, 26, 33, 36, 42], [8, 15, 16, 17, 23, 26, 33, 42], [11, 12, 17, 18, 21, 23, 29, 33, 36, 41, 42, 47], [14, 15, 17, 18, 20, 22, 26, 29, 31, 33, 36, 47]],z,rt,t);
 for(i=1,#sets,z=prod(j=1,#sets[i],Mod(x+sets[i][j],P));
  if(!nfeltissquare(nf,z,&rt),error("unsatisfied relation"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,P)^2!=z,error("identity"));print("SQUARE ",vector(4,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
