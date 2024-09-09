%% FP cont in eps, 
p=spcontini('2/rw1','fpt1',7,'2/fpc1'); huclean(p); plotsol(p); p.sw.bprint=2; 
p.nc.del=1e-5; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=spjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.file.smod=2; p.sw.verb=2; 
p.usrlam=[0.5 0.75]; p.nc.dsmax=0.3; p.sw.bifcheck=0; p=cont(p,10); 
%% HP cont in del=par(7); (increase domain size) 
p=hpcontini('2/rw1','hpt3',7,'2/hpc1'); plotsol(p);
p.nc.del=1e-3; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=[0.5 0.75]; p.file.smod=10; p.nc.dsmax=0.01; 
p=cont(p,130); 
%% HP cont in del=par(7); (increase domain size) 
p=hpcontini('2/rw1','hpt4',7,'2/hpc2'); plotsol(p);
p.nc.del=1e-3; p.nc.tol=1e-6; p.sw.spjac=1; p.sw.jac=1;
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=[0.5 0.75]; p.file.smod=10; p.nc.dsmax=0.01; p=cont(p,80); 
%% plotbra of FPC, HPC
f=5; c=1; mclf(f); 
plotbra('2/fpc1','pt10',f,c,'cl','k'); 
plotbra('2/hpc1',f,c,'cl','r'); 
plotbra('2/hpc2',f,c,'cl','b','lab', [36,70]); 
axis([0.48 1 0 1.25]); xlabel('\epsilon'); grid on; box on; 
%% HP cont exit, with cont in both directions, del=0.75, left HP  
p=hpcontexit('2/hpc1','pt62','2/phpc1'); % puts the HP into dir 'PostHPcont1' 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.05; p.sw.bifcheck=2; 
p=cont(p,25); p=loadp('2/phpc1','pt0','2/phpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% HP cont exit, with cont in both directions, del=0.75, right HP  
p=hpcontexit('2/hpc2','pt36','2/phpc2'); 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.05; p.sw.bifcheck=2; 
p=cont(p,25); p=loadp('2/phpc2','pt0','2/phpc2b'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% HP cont exit, with cont in both directions, del=0.65, right HP  
p=hpcontexit('2/hpc2','pt50','2/phpc2d'); 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.025; p.sw.bifcheck=2; 
p=cont(p,25); p=loadp('2/phpc2d','pt0','2/phpc2db'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% HP cont exit, with cont in both directions, del=0.5, right HP  
p=hpcontexit('2/hpc2','pt70','2/phpc2e'); 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.05; p.sw.bifcheck=2; 
p=cont(p,25); p=loadp('2/phpc2e','pt0','2/phpc2eb'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% branches, del=0.75
f=3; c=9; mclf(f); ax=[-0.2 1.5 0 1.45]; ylab='max(u_1)';
%plotbra('2/phpc1',f,c,'cl','r'); plotbra('2/phpc1b',f,c,'cl','r'); 
plotbra('2/rw1','pt50',f,c,'fp',0,'cl',p2pc('r1')); 
plotbra('2/phpc2',f,c,'cl','b','lsw',0); plotbra('2/phpc2b',f,c,'cl','b','lsw',0);
xlabel('r'); ylabel(ylab); axis(ax); grid on; box on; 
axis([0.02 1.4 0.47 1.2]); 
%% branches, del=0.65, stab range no longer given by [hp1c,hp2c]
% interaction with other modes has happened! 
f=3; c=9; mclf(f); 
plotbra('2/phpc2d',f,c,'cl','b','hplab',1); plotbra('2/phpc2db',f,c,'cl','b'); 
%% branches, del=0.5, stab range (left of hp2c) still exists
f=3; c=9; mclf(f); 
plotbra('2/phpc2e',f,c,'cl','b','labi',0); 
plotbra('2/phpc2eb',f,c,'cl','b','hplab',1); 
xlabel('r'); ylabel(ylab); grid on; box on; 
axis([0.02 0.4 0.25 0.8]); 
%% soln plots 
plotsol('2/hpc2','pt36'); view(vi);nolti; colorbar; pause
plotsol('2/hpc2','pt70'); view(vi);nolti; colorbar; pause 
plotsol('2/hpc1','pt62'); view(vi);nolti; colorbar; pause
plotsol('2/hpc1','pt122'); view(vi);nolti; colorbar; 
%% hoswibra after HPcont (expensive, ca 90s for precon, then ca 20min for 10 steps)
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('phpc2','hpt1',0.04,4,'mrw3cnew',aux); p.sw.verb=2;
p.file.smod=2; p.hopf.flcheck=0; p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; 
p.sw.bifcheck=0; p.nc.dsmax=0.1; p.hopf.ax='unif'; p.file.smod=1; p.plot.bpcmp=9; 
bw=2; beltol=1e-6; belimax=5; % border-width, bel-parameters 
p=setbel(p,bw,beltol,belimax,@lssAMG); p.nc.tol=1e-5; p=cont(p,10);