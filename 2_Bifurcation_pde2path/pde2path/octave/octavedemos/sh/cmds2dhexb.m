keep pphome; close all; p=[]; % localized hexagons for SH on long rectangle 
%% init and zero-branch 
lx=4*pi; nx=round(20*lx); ly=2*pi/sqrt(3); lam=-0.001; nu=1.3; par=[lam; nu];  sw.sym=2; sw.ref=6; nx=3; 
ndim=2; p=shinit(p,nx,lx,ly,ndim,par,sw); p.np, dir='2d_4'; p=setfn(p,dir); p.sol.ds=0.005; 
p.nc.dsmin=0.005; p.sol.dsmax=0.05; p.sw.bifcheck=2; p.np, p=cont(p,5); 
%% qswibra at 0, needs low isotol due to poor grid, then also reset some switches  
aux=[]; aux.m=3; aux.isotol=1e-16; p0=qswibra(dir,'bpt1',aux); 
p0.pm.resfac=1e-4; p0.nc.dsmin=0.03; p0.sw.bifcheck=2; p0.nc.mu2=0.05;
%% hex via seltau 
p=seltau(p0,1,'2dh1b',2); p.sol.ds=-0.05; p=cont(p,20); 
%% gentau for stripes 
p=gentau(p0,[0 1  0],'2ds1b'); p.sol.ds=-0.01; p=cont(p,30); 
%% branch plots
fnr=3; figure(fnr); cmp=3; clf; 
plotbra('2dh1b','pt15',fnr,cmp,'cl','r','lab',20); %,'fancy',1); 
plotbra('2ds1b','pt30',fnr,cmp,'cl','b','lab',20); %,'fancy',1); 
xlabel('\lambda'), ylabel('||u||_*'); axis([-0.22 0.1 0 0.8]); box on; 
%% Zoom 
fnr=3; figure(fnr); cmp=3; clf; 
plotbra('2Dh1b','pt40',fnr,cmp,'cl','r'); 
plotbra('2Ds1b','pt40',fnr,cmp,'cl','b','bplab',16); 
xlabel('\lambda'), ylabel('||u||_*'); axis([-0.05 0.15 0.55 0.65]); box on; 
%% soln plots 
spl('2dh1b','pt15'); spl('2ds1b','pt30'); 