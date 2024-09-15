%% cGL on disk, FPC and HPC in eps 
%% FP cont in eps, 
p=spcontini('rw1','fpt1',7,'fpc1'); p.fuha.outfu=@hpcbra; huclean(p); plotsol(p); 
p.nc.del=1e-5; p.nc.tol=1e-6; 
p.sw.spjac=1; p.sw.jac=1; 
p.sw.jac=0; p.sw.spjac=0; p.nc.njthresh=1e-3; p.nc.njthreshsp=1e3; 
[Ja,Jn]=spjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.file.smod=2; p.sw.verb=2; 
p.usrlam=[0.5 0.75]; p.nc.dsmax=0.3; p.sw.bifcheck=0; p=cont(p,10); 
%% HP cont of HP3 in eps=par(7); (increase domain size) 
p=hpcontini('rw1','hpt3',7,'hpc1'); p.fuha.outfu=@hpcbra; huclean(p); plotsol(p); 
p.nc.del=1e-3; p.nc.tol=1e-6; 
p.sw.spjac=1; p.sw.jac=1; 
%p.sw.jac=1; p.sw.spjac=0; p.nc.njthresh=1e-3; p.nc.njthreshsp=1e6; % needed large! 
%[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=[0.5 0.75]; p.file.smod=10; p.nc.dsmax=0.01; 
p=cont(p,130); 
%% HP cont in del=par(7); (increase domain size) 
p=hpcontini('rw1','hpt4',7,'hpc2'); p.fuha.outfu=@hpcbra; plotsol(p);
p.nc.del=1e-3; p.nc.tol=1e-6; 
p.sw.jac=1;p.sw.spjac=1; 
%p.sw.jac=1; p.sw.spjac=0; p.nc.njthresh=1e-3; p.nc.njthreshsp=1e6; % needed large! 
[Ja,Jn]=hpjaccheck(p); pause 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); p.sw.bprint=p.nc.ilam(2); 
p.usrlam=[0.5 0.75]; p.file.smod=10; p.nc.dsmax=0.01; p=cont(p,75); 
%% plotbra of FPC, HPC, r
f=5; c=1; mclf(f); plotbra('fpc1','pt10',f,c,'cl','k'); 
plotbra('hpc1',f,c,'cl','r'); plotbra('hpc2',f,c,'cl','b','lab', [36,70]); 
axis([0.48 1 0 1.25]); xlabel('\epsilon'); grid on; box on; 
%% plotbra of HPC, omega
f=5; c=10; mclf(f); plotbra('hpc1','pt130',f,c,'cl','r'); 
plotbra('hpc2','pt70',f,c,'cl','b'); 
ylabel('\omega');xlabel('\epsilon'); grid on; box on; 
%% HP cont exit, with cont in both directions, del=0.75, left HP  
p=hpcontexit('hpc1','pt62','phpc1',@hobratw); % puts the HP into dir 'PostHPcont1' 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.05; p.sw.bifcheck=2; 
p=cont(p,25); p=loadp('phpc1','pt0','phpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% HP cont exit, with cont in both directions, del=0.75, right HP  
p=hpcontexit('hpc2','pt36','phpc2',@hobratw); 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.05; p.sw.bifcheck=2; 
p=cont(p,25); p=loadp('phpc2','pt0','phpc2b'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% HP cont exit, with cont in both directions, del=0.65, right HP  
p=hpcontexit('hpc2','pt50','phpc2d',@hobratw); 
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.025; p.sw.bifcheck=2;
p.nc.neig=100; p=cont(p,25); p=loadp('phpc2d','pt0','phpc2db'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% HP cont exit, with cont in both directions, del=0.5, right HP  
p=hpcontexit('hpc2','pt70','phpc2e',@hobratw); p.nc.neig=100;  
mclf(2); p.sw.spcalc=1; p.plot.bpcmp=0; p.nc.dsmax=0.05; p.sw.bifcheck=2; p.file.smod=2; 
p=cont(p,25); p=loadp('phpc2e','pt2','phpc2eb'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% branches, eps=0.75
f=3; c=9; mclf(f); [0.02 1.4 0.47 1.2]; ylab='max(u_1)';
plotbra('rw1','pt50',f,c,'fp',0,'cl',p2pc('r1')); 
plotbra('phpc2',f,c,'cl','b','lsw',0); plotbra('phpc2b',f,c,'cl','b','lsw',0);
xlabel('r'); ylabel(ylab); axis(ax); grid on; box on; 
%% branches, del=0.65
f=3; c=9; mclf(f); plotbra('phpc2d','pt26',f,c,'cl','b','hplab',1); 
plotbra('phpc2db','pt12',f,c,'cl','b'); 
%% branches, del=0.5, 
f=3; c=9; mclf(f); ylab='max(u_1)'; plotbra('phpc2e','pt26',f,c,'cl','b');
plotbra('phpc2eb','pt12',f,c,'cl','b'); plotbra('mrw3e',f,c,'cl',p2pc('v2')); 
plotbra('mrw2e','pt4',f,c,'cl',p2pc('v2')); 
xlabel('r'); ylabel(ylab); grid on; box on; axis([0.02 0.4 0.25 0.8]); 
%% soln plots 
plotsol('hpc2','pt36'); view(vi);nolti; colorbar; pause
plotsol('hpc2','pt70'); view(vi);nolti; colorbar; pause 
plotsol('hpc1','pt62'); view(vi);nolti; colorbar; pause
plotsol('hpc1','pt122'); view(vi);nolti; colorbar; 
%% hoswibra after HPcont (expensive, ca 90s for precon, then ca 20min for 10 steps)
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('phpc2eb','hpt1',0.04,4,'mrw3e',aux); p.sw.verb=2;
p.file.smod=1; p.hopf.flcheck=0; p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; 
p.sw.bifcheck=0; p.nc.dsmax=0.1; p.hopf.ax='unif'; p.plot.bpcmp=9; 
bw=2; beltol=1e-6; belimax=5; % border-width, bel-parameters 
p=setbel(p,bw,beltol,belimax,@lssAMG); pause; p.nc.tol=1e-5; p=cont(p,10);
%% hoswibra after HPcont (expensive, ca 90s for precon, then ca 20min for 10 steps)
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; aux.qfh=@qfh; 
aux.qfhder=@qfhjac; p=hoswibra('phpc2e','hpt1',0.04,4,'mrw2e',aux); p.sw.verb=2; pause 
p.file.smod=1; p.hopf.flcheck=0; p.sw.bifcheck=0; p.hopf.ilam=6; p.nc.ilam=1; 
p.sw.bifcheck=0; p.nc.dsmax=0.1; p.hopf.ax='unif'; p.plot.bpcmp=9; 
bw=2; beltol=1e-6; belimax=5; % border-width, bel-parameters 
p=setbel(p,bw,beltol,belimax,@lssAMG); pause; p.nc.tol=1e-5; p=cont(p,10);