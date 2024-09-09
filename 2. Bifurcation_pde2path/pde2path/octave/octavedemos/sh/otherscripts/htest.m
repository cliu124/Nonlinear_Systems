%% testing 2D-H, init and zero-branch 
lx=2*pi; nx=round(8*lx); ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=1.3; par=[lam; nu];  
p=shinit(p,nx,lx,ly,ndim,par); p=setfn(p,'hh'); p.fuha.outfu=@shbra2d; 
p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.05; 
p.file.smod=1; p.sw.bifcheck=2; p=cont(p,10); 
%% hex det. by QBE, with increased soltol
aux=[]; aux.soltol=1e-12; p0=qswibra('hh','bpt1',aux); p0.nc.dsmin=0.05; 
p0.u(p0.nu+2)=2; % choose different nu for checking 
%%
p=seltau(p0,1,'t1',2); p.sol.ds=-0.05; p.nc.tol=1e-6; p=pmcont(p,40); 
%p=gentau(p0,1,'t2'); p.sol.ds=-0.05; p.nc.dsmin=0.05; p=pmcont(p,20); 
%%
plotbra('t1','pt40',3,6,'cl','r'); plotbra('t1','pt20',3,-1,'cl','r'); 
plotbra('t2','pt20',3,6,'cl','m'); plotbra('t1','pt20',3,-1,'cl','r'); 
%%
%p=loadp('2DH8loc','pt30'); p.nx=126; plotsol(p);
p=loadp('2DH8f','pt80'); p.nx=126; plotsol(p);
%p=loadp('t2','pt10'); plotsol(p); p.nx=50; 
geth2(p,p.u)