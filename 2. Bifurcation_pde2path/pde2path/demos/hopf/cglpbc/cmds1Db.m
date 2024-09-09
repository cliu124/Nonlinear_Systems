%% HP cont and FP cont, 1D, with nq>0, here nq=1 (transl.PC)
% HP cont, par(4)=c3 as new primary 
p=hpcontini('1dtw1b','hpt2',4,'hpc1'); huclean(p); p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=2; 
p.nc.del=1e-3; % slightly increase stepsize for FDs,
p.nc.njthresh=1e-1; p.nc.njthreshsp=1e3; % threshholds for numjac 
p.sw.jac=1;  p.sw.spjac=1; % spjac=0,1 almost equal 
p.sw.jac=0; p.sw.spjac=1;  % uncomment to check numjac 
[Ja,Jn]=hpjaccheck(p); Jd=abs(abs(Ja-Jn)); e1=max(max(Jd)); mclf(10); spy(Jd>0.1*e1); pause
p.nc.dsmax=0.1; p.sol.ds=-0.01; p.nc.tol=1e-6; p.file.smod=2; p.usrlam=[-2,-1.5]; 
p.sw.jac=1; p.sw.spjac=0; p.nc.tol=1e-6; % jac=0 and spjac=2 needs larger tol (for initial step) 
tic; p=cont(p,20); toc
%% branch plot of HPcont 
fn=4;cmp=1; mclf(fn); mclf(fn); plotbra('hpc1','pt10',fn,cmp)
%% HP cont exit 
p=hpcontexit('hpc1','pt17','phpc1'); % puts the HBP into dir 'PostHPcont1' 
mclf(1); p.sol.ds=0.02; p.sw.spcalc=1; p.plot.bpcmp=9; p=cont(p,20);  % cont, in next line in other direction 
p=loadp('phpc1','pt2','phpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% hoswibra after HPcont 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=40; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('phpc1','hpt1',0.04,4,'phc1-1',aux); p.sw.verb=1; p.plot.bpcmp=9;
p.file.smod=2; p.hopf.flcheck=0; p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; p=cont(p,20);
%% FP cont, par(4)=c3 as new primary (change 2nd line to check spjac)  
p=spcontini('1dtw1b','fpt1',4,'fpc1'); huclean(p); p.sw.bprint=2; plotsol(p); 
p.nc.del=1e-3; 
p.nc.njthresh=1e-2; p.nc.njthreshsp=1e4; 
p.sw.jac=1; p.sw.spjac=1; % jac=0 or jac=1 no problem 
p.sw.jac=0; p.sw.spjac=0; % jac=0 or jac=1 no problem 
[Ja,Jn]=spjaccheck(p); Jd=abs(abs(Ja-Jn)); e1=max(max(Jd)); mclf(10); spy(Jd>0.4*e1); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.nc.tol=1e-6; p.usrlam=[-2,-1.5]; 
p.file.smod=2; p.sw.bifcheck=0; p.sw.verb=2; p.nc.dsmax=0.3; tic; p=cont(p,20); toc
%% FP cont exit 
p=spcontexit('fpc1','pt12','pfpc1'); p.plot.bpcmp=0; mclf(2); 
p.sw.spcalc=1; p.sw.bifcheck=2; p=cont(p,15); 
p=loadp('pfpc1','pt2','pfpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,15); % other direction
%% BD plot (compare old and new TW and mTW branches) 
fn=3; cmp=9; mclf(fn); sw=2; 
switch sw; case 1; cmp=9; ym=1.7; ylab='max(u_1)'; case 2; cmp=11; ym=1.1; ylab='||u_1||_2'; end 
plotbra('1dtw1b','pt60',fn,cmp,'fp',0,'cl','k','lsw',0); 
plotbra('1dtw1bs1','pt30',fn,cmp,'fp',0,'cl',p2pc('v1'),'lab',[20]); 
plotbra('phpc1','pt20',fn,cmp,'fp',0,'cl',p2pc('g1'),'lsw',0);
plotbra('phpc1b','pt30',fn,cmp,'fp',0,'cl',p2pc('g2'),'lsw',0,'lp',25);
plotbra('phpc1-1','pt20',fn,cmp,'fp',0,'cl',p2pc('v2'),'lab',[12]); 
axis([0 2 0 ym]); ylabel(ylab); grid on; box on
%% mTW-plots
mclf(1); aux=[]; aux.mper=2; aux.vi=[10,40]; aux.cb=0; 
lframeplot('1dtw1bs1','pt20',11,1,aux); pause 
lframeplot('phpc1-1','pt12',11,1,aux); 
%% TW-plots
mclf(1); plotsol('1dtw1b','fpt1',1,1,1); pause 
plotsol('1dtw1b','hpt2',1,1,1); pause 
plotsol('fpc1','pt16',1,1,1); pause 
plotsol('phpc1','hpt1',1,1,1); pause 
twplot('phpc1','hpt1',1); shading interp; view([10,60]); 
%% branch plot of HPC and FPC
fn=4;cmp=1; mclf(fn); mclf(fn); plotbra('hpc1','pt20',fn,cmp); 
plotbra('fpc1','pt20',fn,cmp,'tyun','--','fms',0,'lab',12); 
axis([-2.1 -1 -0.1 1.3]); 