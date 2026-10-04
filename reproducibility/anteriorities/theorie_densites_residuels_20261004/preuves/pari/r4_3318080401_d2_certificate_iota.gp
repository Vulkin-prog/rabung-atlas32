default(parisizemax,512000000);
default(parisize,128000000);
run()={my(fp=(576)*y^0+(0)*y^1+(49)*y^2+(0)*y^3+(1)*y^4,nf=nfinit(fp),rt,t);if(!nfeltissquare(nf,-1,&rt),error("missing i"));
 t=lift(nfbasistoalg(nf,rt));print("IOTA ",vector(4,j,Str(polcoef(t,j-1))));print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
