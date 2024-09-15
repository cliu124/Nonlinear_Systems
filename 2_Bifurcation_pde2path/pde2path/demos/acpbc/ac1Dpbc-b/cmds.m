%% demo for pBC, 1D with translational invariance 
close all; keep pphome; 
%% init and cont of trivial branch 
p=[]; par=[1 -0.3 1 0]; % here par(4)=s=wave-speed 
p=acinit(p,pi,100,par); p=setfn(p,'tr'); p.nc.dsmax=0.2; p=cont(p,30); 
%% BP1, simple
p=swibra('tr','bpt1','b1',-0.1); p=cont(p); 
%% BP2, double, but just swibra works and simply selects one of the translates. 
p=swibra('tr','bpt2','b2',-0.1); p=cont(p,2); % two cont. steps without PC 
p.nc.nq=1; p.fuha.qf=@qf; p.fuha.qfder=@qfder; % switch on the PC, 
p.nc.ilam=[2 4]; % free the wave speed s as a dummy parameter 
p=cont(p,30); fprintf('s=%g\n',p.u(p.nu+4)); 
%% bifurcation diagram and solution plotting 
f=3; c=0; figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('tr',f,c,'cl','b','lsw',0); plotbra('b1',f,c,'cl','k','lsw',0); 
plotbra('b2',f,c,'cl','r','lab',11); plotsol('b2','pt11'); 
%% 2dim kernel at further BPs, use cswibra to compute and display kernel 
aux=[]; aux.besw=0; p0=cswibra('tr','bpt2',aux); 
%% select 2nd bif direction 
p=gentau(p0,[0 1],'b2alt'); p.sol.ds=-0.1;  p=cont(p,2); % 2 initial steps without PC
p.nc.nq=1; p.fuha.qf=@qfx; p.fuha.qfder=@qfxder; p.nc.ilam=[2 4]; % switch on PC 
p=cont(p,30); fprintf('s=%g\n',p.u(p.nu+4)); plotsol('b2alt','pt9'); % cont and plot