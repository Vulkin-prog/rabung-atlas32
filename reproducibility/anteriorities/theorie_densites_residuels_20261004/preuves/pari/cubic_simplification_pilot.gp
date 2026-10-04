default(parisizemax,512000000);
default(parisize,64000000);
run()={my(P=x^3-197*x^2-9264*x-90252,nf=nfinit(P),z=Mod(x,P));
 print("FACTOR1 ",factor(3506623*x^2+142335035*x+1332041118));
 print("FACTOR2 ",factor(188*x^2+7631*x+71415));
 print("RATIO ",lift((111891453912+11956142940*z+294556332*z^2)/(116257906800+12422657520*z+306048960*z^2)));
 print("END 1");};
iferr(run(),err,print("ERROR ",err);quit(1));quit;
