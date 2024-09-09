%% commands for SH 2nd order sys on perturbed (0,2\pi)^2 square 
keep pphome; close all; p=[]; 
%% init and zero-branch 
lx=2*pi; nx=round(8*lx); ly=lx; ndim=2; lam=-0.001; nu=0.5; par=[lam; nu];  
hmax=0.25; del=0.35; nbd=12; p=shinitpsq(p,nbd,del,hmax,par); p.np 
p.nc.lammax=0.4; p.nc.dsmax=0.05; p.usrlam=0.3; 
p=setfn(p,'2Dp'); huclean(p); p=cont(p,20); 
%% first 2 BPs are now simple  
p=swibra('2Dp','bpt1','mm1'); p=cont(p,20); 
p=swibra('2Dp','bpt2','mm2',-0.05); p=cont(p,20); 
%% stripes by 2ndary bif from 'mm1' 
p=swibra('mm1','bpt1','pm1',-0.1); p=cont(p,20); 
%p=swibra('mm1b','bpt1','pm1b',0.1); p=cont(p,20); 
%% 3rd BP still double, use cswibra 
aux.m=2; p0=cswibra('2Dp','bpt3',aux); 
p0.nc.dsmax=0.1; p0.nc.dsmin=0.05; p0.sol.ds=0.08; 
%% one 'pure' mode and one mixed mode, pmcont important for pure mode
p=seltau(p0,3,'pm3'); p.sol.ds=-0.05; p.pm.mst=10; p.pm.resfac=1e-6; p=pmcont(p,20);
p=seltau(p0,1,'mm3'); p.sol.ds=-0.03; p=cont(p,20); 
%%
p=seltau(p0,3,'pm3b'); p.sol.ds=0.05; p.pm.mst=10; p.pm.resfac=1e-6; p=pmcont(p,20);
p=seltau(p0,1,'mm3b'); p.sol.ds=0.03; p=cont(p,20); 
%% BD 
f=3; c=0; figure(f); clf; plotbra('mm1',f,c,'cl',p2pc('b1')); 
plotbra('mm1b',f,c,'cl',p2pc('b1')); 
plotbra('mm2',f,c,'cl',p2pc('b3')); 
plotbra('pm1',f,c,'cl','k'); 
plotbra('pm1b',f,c,'cl',p2pc('gr1')); 
plotbra('pm3',f,c,'cl',p2pc('r1'));  plotbra('mm3',f,c,'cl',p2pc('r3'));
xlabel('\lambda'); ylabel('||u||_2'); box on; axis([-0.03 0.4 0 3]); 
%% soln plots 
plotsol('mm1','pt16'); nolab; xticks([0 5]); yticks([0 5]); pause; 
plotsol('mm1b','pt13'); nolti; pause
plotsol('mm2','pt12'); nolti; pause
plotsol('pm1','pt11'); nolti; pause; 
plotsol('pm1b','pt7'); nolti; pause; 
%%
plotsol('pm3','pt9'); nolti; pause
plotsol('mm3','pt9'); nolti