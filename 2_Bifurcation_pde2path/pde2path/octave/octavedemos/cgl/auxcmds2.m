%% C1: plotting of precomputed multipliers 
h=plotfloq('1db1','pt39'); figure(7); set(gca,'fontsize',16); 
%% C2: a posteriori compute and plot multipliers
aux=[];aux.wnr=6; aux.nfloq=10; 
[muv1,muv2]=floqap('1db2','pt10',aux); axis tight;
%% C3: this cell only if percomplex has been mexed
aux.wnr=8; [muv1,muv2]=floqpsap('1db2','pt5',aux);
%% C4: time-integration, preparations 
p=loadp('1db1','pt20'); hoplot(p,1,1); dir='stab1d1'; 
p.u(1:p.nu)=p.hopf.y(1:p.nu,1); u0=p.u(1:p.nu); p=setfn(p,dir);
ts=[]; t0=0; npp=50; nt=200; pmod=50; smod=5; tsmod=1; nc=0; 
%% C5: time-integration (repeat if necessary) 
[p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); 
figure(4); clf; plot(ts(1,:), ts(2,:)); % plot values at selected point
figure(5); clf; plot(ts(1,:), ts(3,:)); % plot difference in norm
%% C6: x-t plot; see in ts if there's something interesting after np  
% periods, then plot around there
si=3*npp; incr=25; nt=5*npp/smod; wnr=2; cmp=1; vv=[30,70]; nt=15;
tintplot1d(dir,si,incr,nt,wnr,cmp,vv); 