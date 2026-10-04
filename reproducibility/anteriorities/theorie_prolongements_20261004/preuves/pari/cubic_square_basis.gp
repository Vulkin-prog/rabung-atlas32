default(parisize,128000000);default(parisizemax,512000000);setrand(1);
run()={my(P=x^3+163*x^2+8010*x+121176,nf,sets=[[22, 30, 33, 37], [36, 37, 46, 52], [2, 30, 38, 42], [2, 24, 33, 44], [27, 28, 36, 45, 46, 47], [17, 24, 34, 36, 37, 42], [24, 25, 34, 42, 48, 51], [9, 24, 27, 28, 42, 45], [12, 25, 36, 38, 48, 49], [10, 24, 36, 39, 46, 49], [0, 2, 17, 30, 33, 45], [2, 27, 37, 45, 46, 49]],z,t,rt);nf=nfinit(P);
 for(i=1,#sets,z=prod(j=1,#sets[i],Mod(x+sets[i][j],P));
  if(!nfeltissquare(nf,z,&rt),error("not a square"));t=lift(nfbasistoalg(nf,rt));
  if(Mod(t,P)^2!=z,error("identity"));print("ROOT ",vector(3,j,Str(polcoef(t,j-1)))));
};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
