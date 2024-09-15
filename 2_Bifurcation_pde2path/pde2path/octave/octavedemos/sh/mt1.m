%% commands for SH 2nd order sys on rect. for hexagons
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=2*pi; nx=round(10*lx); ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=0; par=[lam; nu];  
sw.ref=1; sw.sym=0; % ref=#mesh-refs, sym=1: cc-mesh 
nx=10; p=shinit(p,nx,lx,ly,ndim,par,sw);  p.np
 huclean(p); p=setfn(p,'2Dhex0'); plotsol(p,1,1,0); pause 
p.np, p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.05; 
p.file.smod=1; p.sw.bifcheck=2; p=cont(p,10); 