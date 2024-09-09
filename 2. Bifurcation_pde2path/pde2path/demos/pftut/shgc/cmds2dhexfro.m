% localized hexagons for GCSH on long rectangle 
keep pphome; close all; p=[]; global p2pglob;  p2pglob=[]; 
%% init and zero-branch 
lx=8*pi; nx=round(4.5*lx);ly=4*pi/sqrt(3); lam=-0.001; nu=1.3; ga=1; par=[lam; nu; ga];  
ndim=2; sw.sym=1; p=shinit(p,nx,lx,ly,ndim,par,sw); p=setfn(p,'2Dhex8'); p.sol.ds=0.005; 
p.np, pause 
p.nc.dsmin=0.005; p.sol.dsmax=0.05; p.sw.bifcheck=2; p=cont(p,5); 
%% hex from qswibra, need low isotol due to poor grid
aux=[]; aux.isotol=1e-16; aux.m=4; p0=qswibra('2Dhex8','bpt1',aux); 
%%
p=seltau(p0,1,'2Dh8',2); p.pm.resfac=1e-4; p.sol.ds=-0.01; p.sw.eigssol=1; p.sw.verb=2; 
p.nc.dsmin=0.01; p.nc.dsmax=0.05; p.nc.mu2=0.005; p.nc.neig=40; p=pmcont(p,10); 
p.sw.bifcheck=0; % switch off bif-detec for further steps for speed 
p.nc.ntot=1000; p=cont(p,50); 
%% 2ndary bif to hex-front (lower tol and switch off bifcheck for speed) 
p=swibra('2Dh8','bpt1','2DH8f',0.05); p.sw.bifcheck=0; p.nc.tol=1e-6; pause
p.nc.dsmin=0.005; p.nc.dsmax=0.025; p.nc.ntot=1000; p.pm.mst=8; p=pmcont(p,30);
p.nc.ntot=1000; p=cont(p,50); 
%% branch plots 
fnr=3; figure(fnr); cmp=4; clf; 
plotbra('2Dh8','pt30',fnr,cmp,'cl','b','lp',24); 
plotbra('2DH8f','pt80',fnr,cmp,'cl','r','lp',75, 'lab', [20, 50]); 
xlabel('\lambda'), ylabel('||u||'); box on; 
%% soln plots 
d='2DH8f'; pt='pt20'; plotsol(d,pt); xlabel(''); ylabel(''); title(['gchf/' pt]); pause; 
pt='pt50'; plotsol(d,pt); xlabel(''); ylabel(''); title(['gchf/' pt]);