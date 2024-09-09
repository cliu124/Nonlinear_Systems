%% commands for SH 2nd order sys on SQUARE. Hence only pitchforks, use cswibra  
keep pphome; close all; p=[]; 
%% init and zero-branch 
lx=2*pi; nx=round(8*lx); ly=lx; ndim=2; lam=-0.001; nu=0; par=[lam; nu];  
p=shinit(p,nx,lx,ly,ndim,par); p=setfn(p,'2D/0'); huclean(p); p=cont(p,10); 
%% nu=0, hence, sp (spots) supercrit&unstable, st (stripes) supercrit&stable 
p0=cswibra('2D/0','bpt1'); % cswibra, then reset some parameters
p0.nc.dsmax=0.5; p0.nc.dsmin=0.1; p0.nc.lammax=1; p0.sol.ds=0.1; 
%% inspecting the results of cswibra yields tau1=spots, tau3=stripes, hence 
% select these and call cont. 
p=seltau(p0,1,'2D/1-1aSp'); p=pmcont(p,10); 
p=seltau(p0,3,'2D/1-1aSt'); p=pmcont(p,10);
%% nu=nu_1; indeterminate case, cswibra does not find all solns 
p0.u(p0.nu+2)=sqrt(27/38); aux=[]; aux.hasker=1; % subseq. call, reuse kernel 
aux.isotol=1e-1; p0=cswibra(p0,aux); 
%% nu=0.7; a_1>0, a_2<0, a1>-a2, Sp supercrit & stable, St=supercrit & unstable 
p0.u(p0.nu+2)=0.7;  p0=cswibra(p0,aux); p0.sol.ds=0.01;
p=seltau(p0,1,'2D/1-1bSp'); p=cont(p,5); p=seltau(p0,3,'2D/1-1bSt'); p=cont(p,5);
%% nu=0.8; a_1>0, a_2<0, a1<-a2, both unstab
% now we just reset nu and use the old tangents since quadr. predictor makes no difference 
p0.u(p0.nu+2)=0.8; p=seltau(p0,1,'2D/1-1cSp'); p=cont(p,5); p=seltau(p0,3,'2D/1-1cSt'); p=cont(p,5);
%% nu=1; a_1>0, a_2<0, a1<-a2, both unstab
p0.u(p0.nu+2)=1; p=seltau(p0,1,'2D/1-1dSp'); p=cont(p,5); p=seltau(p0,3,'2D/1-1dSt'); p=cont(p,5);
%% plot BD 
fnr=4; cmp=4; figure(fnr); clf; plotbra('2D/0','bpt1',fnr,cmp,'cl','k');
plotbra('2D/1-1aSp','pt10',fnr,cmp,'cl','b','lab',10); 
plotbra('2D/1-1aSt','pt10',fnr,cmp,'cl','r','lab',10); 
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); box on;
%% plot BD, nu=0.7
fnr=4; figure(fnr); clf; plotbra('2D/0','bpt1',fnr,cmp,'cl','k');
plotbra('2D/1-1bSp','pt5',fnr,cmp,'cl','b'); plotbra('2D/1-1bSt','pt5',fnr,cmp,'cl','r'); 
title('\nu=0.7'); xlabel('\lambda'), ylabel('max(|u|)'); box on; 
%% plot BD 
fnr=4; figure(fnr); clf; plotbra('2D/0','bpt1',fnr,cmp,'cl','k');
plotbra('2D/1-1cSp','pt5',fnr,cmp,'cl','b'); plotbra('2D/1-1cSt','pt5',fnr,cmp,'cl','r'); 
title('\nu=0.8'); xlabel('\lambda'), ylabel('max(|u|)'); box on; 
%% plot BD 
fnr=4; figure(fnr); clf; plotbra('2D/0','bpt1',fnr,4,'cl','k');
plotbra('2D/1-1dSp','pt5',fnr,4,'cl','b'); plotbra('2D/1-1dSt','pt5',fnr,4,'cl','r'); 
title('\nu=1'); xlabel('\lambda'), ylabel('max(|u|)'); box on; 
%% soln plots, incl.Fourier 
d='2D/1-1aSp'; pt='pt10'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[3 3],[d '/' pt]); pause
d='2D/1-1aSt'; pt='pt10'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[3 3],[d '/' pt]); pause
d='2D/1-1bSp'; pt='pt5'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[4 4],[d '/' pt]); pause
d='2D/1-1bSt'; pt='pt5'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[4 4],[d '/' pt]); 
%% 2nd BP, Bif. directions not quite obvious from kernel
p0=cswibra('2D/0','bpt2'); p0.nc.dsmax=0.5; p0.nc.dsmin=0.1; p0.nc.lammax=1; p0.sol.ds=0.05; 
p=seltau(p0,1,'2D/1-2a1'); p=pmcont(p,20); p=seltau(p0,3,'2D/1-2a2'); p=pmcont(p,20);
%% 2nd BP, plot BD 
fnr=4; cmp=4; figure(fnr); plotbra('2D/0','bpt1',fnr,cmp,'cl','k');
plotbra('2D/1-2a1','pt15',fnr,cmp,'cl','g','lab',10); 
plotbra('2D/1-2a2','pt20',fnr,cmp,'cl','m','lab',20); 
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); box on;
%% soln plots 
d='2D1-2a1'; pt='pt5'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[5 5],[d '/' pt]); pause
d='2D/1-2a2'; pt='pt20'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[5 5],[d '/' pt]); 
%% other degenerate cases, cswibra finds nonisol. solns 
p0.u(p0.nu+2)=sqrt(81/146); %p0.u(p0.nu+2)=sqrt(27/38); 
del=0.001; p0.u(p0.nu+2)=p0.u(p0.nu+2)+del; 
aux.isotol=1e-8; aux.hasker=1; p0=cswibra(p0,aux); % subseq. call 