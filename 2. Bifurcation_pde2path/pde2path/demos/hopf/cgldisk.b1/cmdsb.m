%% extra commands, FPcont and HPcont (in c3) with nq>0, here nq=1 (rot.PC)
%% FP cont, 
p=spcontini('rw1b','fpt1',4,'fpc1'); huclean(p); plotsol(p); p.sw.bprint=2; 
p.nc.del=1e-5; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=spjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.file.smod=2; p.sw.verb=2; 
p.nc.dsmax=0.3; p.sw.bifcheck=0; p.usrlam=-2; p=cont(p,20); 
%% FP cont exit 
p=spcontexit('fpc1','pt11','pfpc1'); p.plot.bpcmp=0; mclf(2); 
p.sw.bifcheck=2; p=cont(p,30); 
%% other direction
p=loadp('pfpc1','pt2','pfpc1b'); p=resetc(p); p.sol.ds=-p.sol.ds; p=cont(p,16); 
%% HP cont 
p=hpcontini('rw1b','hpt3',4,'hpc1'); huclean(p); plotsol(p);
p.nc.del=1e-5; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=-2; p.file.smod=2; p.nc.del=1e-3; p.nc.dsmax=0.2; p=cont(p,20); 
%% branch plot of HPcont and FPcont 
fn=3;cmp=1; mclf(fn); mclf(fn); plotbra('fpc1','pt20',fn,cmp,'tyun','--')
plotbra('hpc1','pt20',fn,cmp); 
%% HP cont exit, with cont in both directions 
p=hpcontexit('hpc1','pt11','phpc1'); p.file.hcount=2; % puts the HP into dir 'PostHPcont1' 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.tol=1e-6; p.nc.dsmax=0.1; p.sw.bifcheck=2; 
p.nc.p=cont(p,30); p=loadp('phpc1','pt2','phpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% hoswibra after HPcont (expensive, ca 90s for precon, then ca 20min for 10 steps)
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=30; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('phpc1','hpt1',0.04,4,'mrw3new',aux); p.sw.verb=2;
p.file.smod=2; p.hopf.flcheck=0; p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; 
p.sw.bifcheck=0; p.nc.dsmax=0.1; p.hopf.ax='unif'; p.file.smod=1; p.plot.bpcmp=9; 
bw=2; beltol=1e-6; belimax=5; % border-width, bel-parameters 
p=setbel(p,bw,beltol,belimax,@lssAMG); p.nc.tol=1e-5; p=cont(p,10);
%% compare old and new branches 
fn=3; mclf(fn); sw=1; 
switch sw; 
    case 1; cmp=9; ylab='max(u_1)'; ax=[-0.2 1.5 0 1.45]; 
    case 2; cmp=11; ylab='||u_1||_2'; ax=[-0.2 1.5 0 0.45]; 
end 
plotbra('rw1b','pt40',fn,cmp,'fp',0,'cl',p2pc('r1'),'hplab',3); 
plotbra('rw1b3','pt10',fn,cmp,'fp',0,'cl',p2pc('r3'),'lsw',0); 
plotbra('phpc1','pt24',fn,cmp,'cl',p2pc('g1'),'lsw',0); 
plotbra('phpc1b','pt10',fn,cmp,'cl',p2pc('g2'),'lsw',0); 
%plotbra('pfpc1','pt30',fn,cmp,'cl','b','lsw',0); 
plotbra('mrw3new','pt10',fn,cmp,'fp',0,'cl',p2pc('m'),'lsw',0); 
xlabel('r'); ylabel(ylab); axis(ax); grid on; box on; 
%% soln plots: profiles at HPTs 
vi=[0,90]; plotsol('rw1b','hpt3'); view(vi); colorbar; pause; 
plotsol('phpc1','pt0');view(vi); colorbar
%% soln plots, 1/2 period 
hoaux.pind=[1 6 11 16]; hoaux.lay=[1 4]; hoaux.view=[0 90]; hoaux.ztics=[-0.5 0.5]; 
hoplotf('mrw3new','pt20',1,1,hoaux); pause; hoplotf('mrw3new','pt30',1,1,hoaux);
%% Flowers, at pt10: 2nd period approx 5. At pt30: 2nd period approx 3 
%aux.mr=10; aux.pertol=0.025; aux.tol=1e-2; lev=0.25; lfplottip('mrw3new','pt10',10,[lev -lev],aux); pause
aux.mr=12; aux.pertol=0.025; aux.tol=1e-2; lev=0.25; lfplottip('mrw3new','pt30',10,[lev -lev],aux); 
%% FP cont exit 
p=spcontexit('fpc1','pt14','pfpc1'); p.plot.bpcmp=0; mclf(2); p=cont(p,25); 
%% other direction
p=loadp('pfpc1','pt2','pfpc1b'); p=resetc(p); p.sol.ds=-p.sol.ds; p=cont(p,16); 
%% compare old and new branches, including HP-cont branches 
fn=3; mclf(fn); sw=1; 
switch sw; 
    case 1; cmp=9; ylab='max(u_1)'; 
    case 2; cmp=11; ylab='||u_1||_2'; 
end 
plotbra('rw1b','pt40',fn,cmp,'fp',0,'cl',p2pc('r1'),'lab',30); 
plotbra('phpc1','pt26',fn,cmp,'cl',p2pc('g1')); 
plotbra('phpc1b','pt10',fn,cmp,'cl',p2pc('g2')); 
plotbra('pfpc1','pt22',fn,cmp,'cl',p2pc('v1')); 
plotbra('pfpc1b','pt16',fn,cmp,'cl',p2pc('v2')); 
xlabel('r'); ylabel(ylab); 