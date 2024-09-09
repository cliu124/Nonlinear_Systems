%% commands for SH 2nd order sys on CUBE 
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=pi; ly=lx; lz=lx; ndim=3; lam=-0.001; nu=0; par=[lam; nu];  
% choose discretization, including type; 
% nx=3; sw.sym=1; sw.ref=1;
nx=6; sw.sym=1; sw.ref=1; % sym=1: cc-cube, ref=# of (uniform) refinement steps 
%nx=10; sw.sym=0; sw.ref=1; % alternative 
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p=setfn(p,'3D0'); plotsol(p,1,1,3); p.np, pause 
p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; %p.fuha.lss=@lssAMG; 
p.sw.bifcheck=2; p.pm.resfac=1e-3; p=cont(p,5); 
%% cswibra: low isotol cause else solns identified as isolated 
% returns 8 bif. directions; 
aux=[]; aux.m=3; aux.ds=0.02; aux.pstyle=2; %aux.isotol=1e-16; 
p0=cswibra('3D0','bpt1',aux); 
% reset controls for cont of nontr. branches 
p0.nc.dsmin=0.02; p0.nc.dsmax=0.03; p0.sw.bifcheck=0; p0.pm.resfac=1e-4; p0.pm.mst=6; 
%% select 3 bif. directions, lamellas, tubes, and rhombs
p=seltau(p0,1,'3D0-l'); pause, p=pmcont(p,20); % lamella
p=seltau(p0,4,'3D0-t'); pause, p=pmcont(p,20); % tube
p=seltau(p0,11,'3D0-r'); pause, p=pmcont(p,20); % rhomb
%% plot BD and solns 
fnr=3; cmp=4; figure(fnr); clf; plotbra('3D0-t','pt20',fnr,cmp,'cl','b','lab',10);
plotbra('3D0-l','pt20',fnr,cmp,'cl','r','lab',15); plotbra('3D0-r','pt20',fnr,cmp,'cl','m','lab',15);
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); box on; 
plotsol('3D0-t','pt20'); pause; plotsol('3D0-l','pt20'); pause; plotsol('3D0-r','pt20'); 
%% nu=1.2, all branches subcritical 
p0.u(p0.nu+2)=1.2; p0.nc.dsmax=0.5; p0.nc.dsmin=0.05; p0.sol.ds=0.05; p0.file.smod=2; p0.nc.lammax=0.2; 
p=seltau(p0,1,'3D0c-l'); p=pmcont(p,15); % lamella
p=seltau(p0,4,'3D0c-t'); p=pmcont(p,15); % tube
p=seltau(p0,11,'3D0c-r'); p.nc.tol=1e-6; p=pmcont(p,5); p.nc.tol=1e-8; p=pmcont(p,15); % rhomb 
%% BD and solution plot 
fnr=3; cmp=4; figure(fnr); clf; 
plotbra('3D0c-l','pt15',fnr,cmp,'cl','r');
plotbra('3D0c-t','pt10',fnr,cmp,'cl','b'); 
plotbra('3D0c-r','pt10',fnr,cmp,'cl','m'); 
title('\nu=1.2'); xlabel('\lambda'), ylabel('max(|u|)'); axis([-0.08 0.15 0 1.5]); box on; 
plotsol('3D0c-r','pt10',1,1,2); view(-15,10); 
%% nu=0.75: rhombs supercrit and stable at bif, but continue as tubes! 
huclean(p); p0.u(p0.nu+2)=0.75; 
%p=seltau(p0,1,'3D0b-l');  p=pmcont(p,10); % lamella
%p=seltau(p0,5,'3D0b-t'); p=pmcont(p,10); % tube 
p=seltau(p0,7,'3D0b-r'); pause, p.nc.tol=1e-8; p.pm.resfac=1e-6; p.file.smod=5; p=pmcont(p,10); % rhomb 
%% nu=0.75 plots
fnr=3; cmp=4; figure(fnr); clf; plotbra('3D0b-r','pt5',fnr,cmp,'cl','m');
plotbra('3D0b-t','pt10',fnr,cmp,'cl','b'); plotbra('3D0b-l','pt10',fnr,cmp,'cl','r');
title('\nu=0.75'); xlabel('\lambda'), ylabel('max(|u|)'); box on; 
plotsol('3D0b-r','pt5',1,1,2); 