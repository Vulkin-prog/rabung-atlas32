default(parisizemax,1024000000);
default(parisize,128000000);
setrand(1);
run()={my(P=(1)*x^3+(-197)*x^2+(-9264)*x^1+(-90252)*x^0,fp=(-566104520160)*x^0+(0)*x^1+(4435693201)*x^2+(0)*x^3+(-133202)*x^4+(0)*x^5+(1)*x^6,C,raw,nf,fac,roots);
 
 fp=subst(fp,x,y);nf=nfinit(fp);print("STAGE nfinit done");fac=nffactor(nf,P);
 if(matsize(fac)[1]!=3,error("not split"));roots=vector(3,i,-polcoef(fac[i,1],0));
 print("FIELD ",vector(poldegree(fp)+1,i,Str(polcoef(fp,i-1))));
 print("ROOTS ",vector(3,i,vector(poldegree(fp),j,Str(polcoef(lift(roots[i]),j-1)))));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
