close all; keep pphome; 
%% commands for SH on BCC 
p=[]; lx=sqrt(2)*pi; ly=lx; lz=lx; ndim=3; lam=-0.001; nu=1; par=[lam; nu];  
%nx=2; sw.sym=1; sw.ref=3;  dir='BCC0'; % np=3221 nx=2 => best symmetry 
nx=30; sw.sym=1; sw.ref=2; dir='BCC1'; % np=10000, seems best 
%nx=6; sw.sym=2; sw.ref=1; dir='BCC2'; % np=6006, good, but tubes drift
%nx=2; sw.sym=2; sw.ref=4; dir='BCC3'; % np=24.000 slow 
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p.pdeo.grid.plotFaces; p.np, pause 
p.plot.pstyle=3; 
p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; p=setfn(p,dir); 
p.sw.bifcheck=2; p.pm.resfac=1e-3; p=cont(p,10); 
%% qswibra; 3 tubes as kernel; afterwards set some switches 
aux=[]; aux.m=3; p0=qswibra(dir,'bpt1',aux); 
p0.nc.dsmin=0.05; p0.pm.mst=4;  p0.sw.bifcheck=0; p0.pm.resfac=1e-4; 
p0.file.smod=5; p0.nc.tol=1e-6; p0.sw.foldcheck=0; p0.sw.spcalc=1; 
%% select BCC and cont in both directions. Cold BCC tend to loose symmetry, 
%p=seltau(p0,1,'BCCa',2); p.sol.ds=0.02; p=pmcont(p,15); % hence fewer steps for cold 
p=seltau(p0,1,'BCCb',2); p.sol.ds=-0.02; p=cont(p,40); % hot BCC work well  
%% plot BD 
fnr=3; cmp=3; figure(fnr); clf; plotbra('BCC1',fnr,cmp,'cl','k','lp',3); 
plotbra('BCCa','pt15',fnr,cmp,'cl',p2pc('r1'),'lab',15);
plotbra('BCCb','pt20',fnr,cmp,'cl',p2pc('r2'),'lab',20); 
xlabel('\lambda'), ylabel('||u||'); box on; 
%% solns plots 
v=[-50,25]; 
plotsol('BCCa','pt5'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('BCCb','pt20',1,1,2); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 