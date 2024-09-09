close all; keep pphome; 
%% commands for Schnakenberg on a SC cube 
p=[]; kc=sqrt(sqrt(2)-1); lx=pi/kc; ly=lx; lz=lx; par=[3.26, 0, 60];
nx=40; sw.sym=1; sw.ref=1; dir='3DSC'; % np=6006, still OK
p=schnakinit(p,[lx,ly,lz],nx,par,sw); p.pdeo.grid.plotFaces; p.np, pause 
p.pm.resfac=1e-4; p.sol.ds=-0.1; p=setfn(p,dir); p=cont(p,3); p.sw.bifcheck=0; p=cont(p,10); 
%% cswibra; eigenvectors are rather nice lamellas. Thus, running with low soltol, for 
% lamellas and tubes cswibra simply gives back al=(0 0 1) and al=(0 1 1) as solutions; 
% for rhombs: small corrections. 
aux=[]; aux.soltol=1e-12; aux.del=1e-4; 
aux.besw=0; % only kernel, and only lamellas below (saving time!) 
p0=cswibra('3DSC','bpt1',aux); p0.nc.dsmin=0.1; p0.sw.bifcheck=0; 
p0.pm.resfac=1e-6; p0.pm.mst=10; 
p0.file.smod=5; p0.nc.tol=1e-8; p0.nc.lammin=3; p0.sol.ds=0.01; 
%% lamella's and rhombs work with cont without loosing symmetry, tubes need pmcont, 
% but still drift away near lambda=2.5
%p=seltau(p0,1,'3DSC-l'); pause, p=cont(p,20); % lamellas 
%p=seltau(p0,10,'3DSC-r'); pause, p.nc.tol=1e-5; p.sol.ds=-0.05; p=pmcont(p,20); 
%p=seltau(p0,2,'3DSC-t'); p.nc.tol=1e-6; pause, p.sol.ds=0.05; p=pmcont(p,20); % tubes need pmcont 
%% alternative: use gentau 
p=gentau(p0,[0 0 1],'3D0-l'); p.plot.pstyle=3; p=cont(p,20); % lamella
% p=gentau(p0,[1 -1],'3D0-t'); p=pmcont(p,20); % tube  
%p=gentau(p0,[1 -1 1],'3D0-r'); p=cont(p,20); % rhomb
%% plot BD 
fnr=3; cmp=4; figure(fnr); clf; plotbra('3DSC','bpt1',fnr,cmp,'cl','k','lp',6);
plotbra('3DSC-t','pt20',fnr,cmp,'cl','b');pause
plotbra('3DSC-l','pt20',fnr,cmp,'cl','r'); 
plotbra('3DSC-r','pt20',fnr,cmp,'cl','m','lab',[5 10 20]);
xlabel('\lambda'), ylabel('max(|u|)'); box on; 
%% plot solns
v=[-20,10]; plotsol('3DSC-r','pt5'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('3DSC-t','pt10'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
