default(parisize,128000000);default(parisizemax,512000000);setrand(1);
run()={my(P=x^4+80*x^3+2281*x^2+27240*x+116640,f,fp,bn,fac,roots);
 gettime();f=polcompositum(P,x^2-265,1)[1][1];fp=subst(polredbest(f),x,y);
 print("FIELD ",fp);bn=bnfinit(fp,1);print("BNF_MILLISECONDS ",gettime());
 fac=nffactor(bn,P);roots=vector(4,i,-polcoef(fac[i,1],0));
 print("ROOTS ",roots);print("DISCRIMINANT ",nfdisc(fp));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
