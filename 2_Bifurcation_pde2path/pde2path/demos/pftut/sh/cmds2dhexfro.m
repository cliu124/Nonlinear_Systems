keep pphome; close all; p=[]; % localized hexagons for SH on long rectangle 
%% init and zero-branch 
lx=8*pi; nx=round(5*lx);ly=4*pi/sqrt(3); lam=-0.001; nu=1.3; par=[lam; nu];  
ndim=2; p=shinit(p,nx,lx,ly,ndim,par); p=setfn(p,'2D/hex8'); p.sol.ds=0.005; 
p.nc.dsmin=0.005; p.sol.dsmax=0.05; p.sw.bifcheck=2; p=cont(p,5); 
%% hex from qswibra 
p0=qswibra('2D/hex8','bpt1'); p0.sw.bifcheck=1; 
%%
p=seltau(p0,1,'2D/h8',2); p.pm.resfac=1e-4; p.sol.ds=-0.005; p.nc.dsmin=0.01; 
p.u(p.nu+2)=1.3; p.nc.mu2=0.01; p.nc.neig=40; p=pmcont(p,10); 
p.sw.bifcheck=0; % switch off bif-detec for further steps for speed 
tic; p=pmcont(p,50); toc
%% 2ndary bif to hex-front (lower tol and switch off bifcheck for speed) 
p=swibra('2D/h8','bpt1','2D/H8f',-0.001); p.sw.bifcheck=0; p.nc.tol=1e-6; pause
p=cont(p,2); 
p.nc.dsmin=0.001; p.nc.dsmax=0.05; p.nc.ntot=1000; p.pm.mst=8; p=pmcont(p,120);
%% 2ndary bif to hex-loc 
p=swibra('2D/h8','bpt2','2D/H8loc',0.05); p.sw.bifcheck=0; p.nc.tol=1e-6; 
p.nc.dsmin=0.001; p.nc.dsmax=0.05; p.nc.ntot=1000; p.pm.mst=6; p=pmcont(p,120);
%% branch plots 
fnr=3; figure(fnr); cmp=3; clf; 
plotbra('2D/h8','pt50',fnr,cmp,'cl','b','lp',45); 
plotbra('2D/H8f','pt100',fnr,cmp,'cl','r','lp',150, 'lab', [40, 70]); 
plotbra('2D/H8loc','pt60',fnr,cmp,'cl','m','lp',150, 'lab', [30,40]); 
xlabel('\lambda'), ylabel('||u||_*'); box on; 
%% soln plots 
d='2D/H8f'; pt='pt40'; plotsol(d,pt); xlabel(''); ylabel(''); title(['hf/' pt]); pause; 
pt='pt70'; plotsol(d,pt); xlabel(''); ylabel(''); title(['hf/' pt]);
%%
d='2D/H8loc'; pt='pt30'; plotsol(d,pt); xlabel(''); ylabel(''); title(['hloc/' pt]); pause; 
pt='pt40'; plotsol(d,pt); xlabel(''); ylabel(''); title(['hloc/' pt]);