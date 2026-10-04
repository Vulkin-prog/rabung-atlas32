default(parisizemax,512000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(1)*x^4+(100)*x^3+(3258)*x^2+(34844)*x^1+(26309)*x^0,nf=nfinit(P),sets=[[1, 21, 22], [3, 21, 22, 23, 29, 38, 47, 49], [5, 21, 22, 29, 37, 38, 39, 45, 47], [7, 21, 22, 23, 29, 31, 33, 37, 38, 47], [8, 21, 23, 29, 31, 37, 38], [13, 20, 23, 31, 33, 45], [15, 20, 22, 23, 29, 33, 37, 38, 41, 45, 47], [17, 20, 22, 23, 29, 31, 33, 38, 41, 45]],z,rt,t);
 for(i=1,#sets,z=prod(j=1,#sets[i],Mod(x+sets[i][j],P));
  if(!nfeltissquare(nf,z,&rt),error("unsatisfied relation"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,P)^2!=z,error("identity"));print("SQUARE ",vector(4,j,Str(polcoef(t,j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
