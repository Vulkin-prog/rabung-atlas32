run()={my(P=x^4+80*x^3+2281*x^2+27240*x+116640,nf=nfinit(P),I=[0, 5, 12, 24, 27],z,rt,t);z=prod(j=1,#I,Mod(x+I[j],P));if(!nfeltissquare(nf,z,&rt),error("not square"));t=lift(nfbasistoalg(nf,rt));print("SQUARE ",vector(4,j,Str(polcoef(t,j-1))));print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
