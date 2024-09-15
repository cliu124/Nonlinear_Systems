close all; keep pphome; p=[]; 
%% CH 2D, with FPC and HPC; eps=0.1, init and cont of homogen. branch 
p=[]; m=-1; eps=1/10; lam=0; par=[m eps lam]; lx=[0.5 0.5]; nx=[40 40]; 
p=chinit(p,lx,nx,par); p=setfn(p,'2D'); p.nc.dsmax=0.04; p=cont(p,100); 
%% 1st BP, double, hence cswibra 
p0=cswibra('2D','bpt1',aux); p0.nc.lammax=0.3; 
p=seltau(p0,3,'2D1-sp',3); p.sol.ds=0.1; p=cont(p,50); % spots
p=seltau(p0,1,'2D1-st',3); p.sol.ds=-0.1; p=cont(p,60); % stripes
%% 2ndary bif
p=swibra('2D1-st','bpt4','du2'); p=cont(p,40);  
%% BD plot E over m
figure(3); clf; plotbra('2D','pt30','lsw',0); 
plotbra('2D1-sp','pt50','cl','b','lsw',0); 
plotbra('2D1-st','pt60','cl','r','lsw',0); 
plotbra('du2','pt40','cl','m','lsw',0); 
axis([-0.9 -0.3 0.4 2.2]); ylabel('E_\epsilon'); 
%% sol plots
mypsol('2D1-sp','fpt1'); mypsol('2D1-st','bpt1'); mypsol('du1','pt10'); 
%% continuation in eps of FP1 on 2D1-st: 
figure(2); clf; p=spcontini('2D1-st','fpt1',2,'fp2c'); p.sol.ds=-0.001; plotsol(p); 
p.sw.bifcheck=0; p.nc.dsmax=0.05; p.file.smod=5; p.usrlam=[0.075 0.05]; 
p.nc.del=1e-3; p.nc.tol=1e-6; huclean(p); p=cont(p,20); 
%% BD for fold-cont 
figure(3); clf; plotbra('fp2c','pt20',3,1); xlabel('\epsilon'); 
%% soln plots 
mypsol('fp2c','pt0'); mypsol('fp2c','pt12');
%% FP cont exit 
huclean(p); p=spcontexit('fp2c','pt12','st2a'); p.sw.spcalc=1; p.usrlam=0; 
p.nc.dsmax=0.2; p.nc.lammax=1; p.nc.lammin=-1; p.plot.bpcmp=5;  p=cont(p,30); 
%% continuation in eps of BP4 on 2D1-st: 
figure(2); clf; p=bpcontini('2D1-st','bpt4',2,'bp2c'); p.sol.ds=-0.001; plotsol(p); 
p.sw.bifcheck=0; p.nc.dsmax=0.05; p.file.smod=5; p.usrlam=[0.075 0.05]; 
p.nc.del=1e-2; p.nc.tol=1e-6; huclean(p); p=cont(p,20); 
%% BP cont exit 
huclean(p); p=bpcontexit('bp2c','pt13','pbpc2'); p.sw.spcalc=1; p.usrlam=0; 
p.nc.dsmax=0.2; p.nc.lammax=1; p.nc.lammin=-1; p.plot.bpcmp=5;  p=cont(p,20); 
%% also other direction 
p=loadp('pbpc2','bpt1','pbpc2b'); p.sol.ds=-10*p.sol.ds; p=cont(p,20); 
%% check that BP is fine, i.e., that branch-switching works
p=swibra('pbpc2','bpt1','du2b',0.1); p.nc.lammin=-1; p.plot.bpcmp=5; p=cont(p,10); 
%% plot new branches, where b1 and c1 should be identical, which they are
f=4; c=5; mclf(f); plotbra('2D1-st','pt40',f,c,'cl','r','lsw',0); % old, with eps=0.1; 
plotbra('du2','pt20',f,c,'cl','m','lsw',0); 
plotbra('pbpc2','pt20',f,c,'cl',p2pc('g3'),'lsw',0); % new, with eps=0.075
plotbra('pbpc2b','pt20',f,c,'cl',p2pc('g3'),'lsw',0); 
plotbra('du2b','pt10',f,c,'cl',p2pc('g1'),'lsw',0); 
ylabel('E_\epsilon'); %axis([-0.95 -0.4 0.93 1.1]); 