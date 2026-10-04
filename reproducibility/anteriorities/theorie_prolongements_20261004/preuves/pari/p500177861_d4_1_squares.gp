default(parisizemax,512000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(1)*x^4+(115)*x^3+(4313)*x^2+(57861)*x^1+(224350)*x^0,nf=nfinit(P),sets=[[4, 5, 14, 15, 18, 19, 21, 31, 33, 34, 39, 41, 49, 50], [7, 13, 14, 15, 19, 21, 34, 41, 43, 45], [8, 13, 21, 25, 28, 29, 43, 46, 49, 50], [11, 13, 25, 41, 43, 50]],z,rt,t);
 for(i=1,#sets,z=prod(j=1,#sets[i],Mod(x+sets[i][j],P));
  if(!nfeltissquare(nf,z,&rt),error("unsatisfied relation"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,P)^2!=z,error("identity"));print("SQUARE ",vector(4,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
