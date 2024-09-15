%% cGL on disk, FPC and HPC in c3 
%% FP cont, 
p=spcontini('rw1','fpt1',4,'fpc1a'); mclf(2); plotsol(p); p.sw.bprint=2; 
p.nc.del=1e-5; p.nc.tol=1e-6; 
p.sw.jac=1; p.sw.spjac=1; 
%p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e-3; p.nc.njthreshsp=1e3; 
%[Ja,Jn]=spjaccheck(p); pause % uncomment to check impl.of spjac.m 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.file.smod=2; p.sw.verb=2; 
p.nc.dsmax=0.3; p.sw.bifcheck=0; p.usrlam=-2; p=cont(p,20); 
%% FP cont exit 
p=spcontexit('fpc1a','pt11','pfpc1a'); p.plot.bpcmp=0; mclf(2); 
p.sw.bifcheck=2; p.nc.dsmax=0.1; p=cont(p,40); 
p=loadp('pfpc1a','pt2','pfpc1ab'); p=resetc(p); p.sol.ds=-p.sol.ds; % other direction
p=cont(p,16); 
%% HP cont of HP3
p=hpcontini('rw1','hpt3',4,'hpc1a'); plotsol(p); p.fuha.outfu=@hpcbra; 
p.nc.del=1e-2; p.nc.tol=1e-6; 
p.sw.spjac=1; p.sw.jac=1;
p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e-3; p.nc.njthreshsp=1e6; % needed large! 
[Ja,Jn]=hpjaccheck(p); pause  % uncomment to check impl.of spjac.m 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=-2; p.file.smod=10; p.nc.dsmax=0.2; p=cont(p,120); 
%% HP cont of HP4 
p=hpcontini('rw1','hpt4',4,'hpc2a'); plotsol(p);p.fuha.outfu=@hpcbra;
p.nc.del=1e-3; p.nc.tol=1e-6; 
p.sw.spjac=1; p.sw.jac=1;
%p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e-3; p.nc.njthreshsp=1e6; % needed large! 
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=-2; p.file.smod=2; p.nc.dsmax=0.2; p=cont(p,15); 
%% branch plot of HPcont and FPcont 
f=4;cmp=1; mclf(f); plotbra('fpc1a','pt16',f,cmp,'tyun','--')
plotbra('hpc1a','pt90',f,cmp,'cl','r'); plotbra('hpc2a','pt16',f,cmp,'cl','b'); 
axis([-2.5 -1 -0.4 1.3]); box on; grid on; 
%% plotbra of HPC, omega
f=5; cmp=10; mclf(f); plotbra('hpc1a','pt90',f,cmp,'cl','r'); 
plotbra('hpc2a','pt16',f,cmp,'cl','b'); 
ylabel('\omega');xlabel('\epsilon'); grid on; box on; 
%% HP cont exit, with cont in both directions 
p=hpcontexit('hpc2a','pt11','phpc2a'); p.file.hcount=2; % puts the HP into dir 'PostHPcont1' 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.tol=1e-6; p.nc.dsmax=0.1; p.sw.bifcheck=2; 
p.nc.p=cont(p,30); p=loadp('phpc2a','pt2','phpc2ab'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
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
plotbra('rw1','pt50',fn,cmp,'fp',0,'cl',p2pc('r1'),'hplab',[3,4]); 
plotbra('mrw2','pt10',fn,cmp,'fp',0,'cl','r','lsw',0); 
plotbra('mrw3','pt10',fn,cmp,'fp',0,'cl','b','lsw',0); 
plotbra('pfpc1','pt40',fn,cmp,'fp',0,'cl',p2pc('g1')); 
plotbra('pfpc1b','pt16',fn,cmp,'fp',0,'cl',p2pc('g1')); 
xlabel('r'); ylabel(ylab); axis(ax); grid on; box on; 
%% soln plots: profiles at HPTs 
vi=[0,90]; plotsol('rw1','hpt3'); view(vi); colorbar; nolti; pause; 
plotsol('rw1','hpt4'); view(vi); colorbar;nolti; pause; 
plotsol('hpc1a','pt52'); view(vi);nolti; colorbar; pause
plotsol('hpc2a','pt11'); view(vi);nolti; colorbar; 