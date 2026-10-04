default(parisizemax,512000000);
default(parisize,128000000);
run()={my(fp=(-566104520160)*y^0+(0)*y^1+(4435693201)*y^2+(0)*y^3+(-133202)*y^4+(0)*y^5+(1)*y^6,nf=nfinit(fp),alpha=Mod((4078473223/17076311)*y^0+(0)*y^1+(-333005/34152622)*y^2+(0)*y^3+(3/34152622)*y^4+(0)*y^5,fp),cs=[2, 3, 5, 11, 17, 19, 59, 763],sets=[[3, 12, 15, 20, 22, 24, 26, 29, 43], [2, 3, 15, 17, 20, 24, 26, 39, 43], [2, 3, 8, 12, 15, 17, 18, 24, 26, 27, 29, 30, 37, 39, 43], [3, 9, 17, 20, 22, 24, 30, 43], [3, 8, 15, 18, 20, 24, 27, 29, 30, 37, 43], [8, 12, 15, 20, 26, 29, 30, 37, 43], [3, 7, 9, 13, 19, 20, 26, 30, 31], [0, 3, 9, 10, 11, 12, 18, 22, 23, 24, 29, 37, 42, 43]],z,rt,t);
 for(i=1,#cs,z=cs[i]*prod(j=1,#sets[i],alpha+sets[i][j]);if(!nfeltissquare(nf,z,&rt),error("spurious common class"));
  t=lift(nfbasistoalg(nf,rt));print("SQUARE ",vector(6,j,Str(polcoef(t,j-1)))));print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
