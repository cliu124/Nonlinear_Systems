%% extra commands, FPcont and HPcont (in c3) with nq>0, here nq=1 (rot.PC)
%% FP cont, 
p=spcontini('2/rw1','fpt1',4,'2/fpc1a'); huclean(p); plotsol(p); p.sw.bprint=2; 
p.nc.del=1e-5; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=spjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.file.smod=2; p.sw.verb=2; 
p.nc.dsmax=0.3; p.sw.bifcheck=0; p.usrlam=-2; p=cont(p,20); 
%% FP cont exit 
p=spcontexit('2/fpc1a','pt11','2/pfpc1a'); p.plot.bpcmp=0; mclf(2); 
p.sw.bifcheck=2; p.nc.dsmax=0.1; p=cont(p,40); 
%% other direction
p=loadp('2/pfpc1a','pt2','2/pfpc1ab'); p=resetc(p); 
p.sol.ds=-p.sol.ds; p=cont(p,16); 
%% HP cont 
p=hpcontini('2/rw1','hpt3',4,'2/hpc1a'); plotsol(p);
p.nc.del=1e-3; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=-2; p.file.smod=10; p.nc.del=1e-3; p.nc.dsmax=0.2; p=cont(p,20); 
%% HP cont 
p=hpcontini('2/rw1','hpt4',4,'2/hpc2a'); plotsol(p);
p.nc.del=1e-3; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=-2; p.file.smod=2; p.nc.del=1e-3; p.nc.dsmax=0.2; p=cont(p,15); 
%% branch plot of HPcont and FPcont 
fn=4;cmp=1; mclf(fn); mclf(fn); plotbra('2/fpc1','pt16',fn,cmp,'tyun','--')
plotbra('2/hpc1a','pt130',fn,cmp,'cl','r'); 
plotbra('2/hpc2a','pt16',fn,cmp,'cl','b'); 
axis([-2.5 -1 -0.4 1.3]); box on; grid on; 
%% HP cont exit, with cont in both directions 
p=hpcontexit('2/hpc2a','pt11','2/phpc2a'); p.file.hcount=2; % puts the HP into dir 'PostHPcont1' 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.tol=1e-6; p.nc.dsmax=0.1; p.sw.bifcheck=2; 
p.nc.p=cont(p,30); p=loadp('2/phpc2a','pt2','2/phpc2ab'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
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
plotbra('2/rw1','pt50',fn,cmp,'fp',0,'cl',p2pc('r1'),'hplab',[3,4]); 
plotbra('2/mrw2','pt10',fn,cmp,'fp',0,'cl','r','lsw',0); 
plotbra('2/mrw3','pt10',fn,cmp,'fp',0,'cl','b','lsw',0); 
plotbra('2/pfpc1a','pt40',fn,cmp,'fp',0,'cl',p2pc('g1')); 
plotbra('2/pfpc1ab','pt16',fn,cmp,'fp',0,'cl',p2pc('g1')); 
xlabel('r'); ylabel(ylab); axis(ax); grid on; box on; 
%% soln plots: profiles at HPTs 
vi=[0,90]; plotsol('2/rw1','hpt3'); view(vi); colorbar; nolti; pause; 
plotsol('2/rw1','hpt4'); view(vi); colorbar;nolti; pause; 
plotsol('2/hpc1a','pt89'); view(vi);nolti; colorbar; pause
plotsol('2/hpc2a','pt11'); view(vi);nolti; colorbar; 