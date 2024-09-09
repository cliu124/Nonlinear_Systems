close all; keep pphome; 
%% Schnakenberg on larger 2D domain with hex lattice (snaking of str2hex front) 
p=[]; kc=sqrt(sqrt(2)-1); lx=6*pi/kc; ly=lx/sqrt(3)/2; par=[3.21, 0, 60]; 
nx=100; ny=round(nx*ly/lx); Lfn='Ll'; % filename for L 
p=schnakinit(p,lx,ly,nx,ny,par,Lfn); % init 
p.pm.resfac=1e-4; p.sol.ds=-0.01; p=setfn(p,'homl');plotsol(p); p.sw.verb=2; p.sw.bifcheck=2; 
p.nc.bisecmax=5; p.np, p=findbif(p,1);
%% hex via qswibra, continue in both directions
p0=qswibra('homl','bpt1'); p0.nc.dsmin=0.1;  %p0.sw.bifcheck=0; pause
p0.nc.bisecmax=4; p0.nc.mu2=0.02;  p0.nc.dsmax=0.4; 
p0=setbelilup(p0,1,1e-4,5,1e-6,200); % ilupack not efficient here 
%% hex+
p=seltau(p0,2,'h+l',2); p.sol.ds=-0.1; p.fuha.innerlss=@lss; % ilupack not efficient here 
p=pmcont(p,20);
%% s+ via gentau
p=gentau(p0,1,'s+l'); p=cont(p,15); 
%% b+ 
p=swibra('s+l','bpt3','b+l',-0.05); p.nc.dsmin=0.05; p.nc.dsmax=0.1; 
p=findbif(p,2); p=cont(p,15);
%% snake of str2hex front
ds=0.02; p=swibra('b+l','bpt1','sn1',-ds); pause; p.nc.dsmin=ds; p.nc.dsmax=2*ds;  
p.sw.bifcheck=0; p=pmcont(p,20);
%% plot BD 
fnr=3; figure(fnr); clf; plotbra('homl',fnr,4,'cl','k','lsw',0);
plotbra('s+l',fnr,4,'cl','b'); 
plotbra('h+l',fnr,4, 'cl','m','lab',20);  
plotbra('b+l',fnr,4, 'cl','r','lab',20);  
plotbra('sn1',fnr,4, 'cl',p2pc('o1'),'lab',20);  
xlabel('\lambda'); ylabel('max/min(u_1)');  
%% plot sol. incl. Fourier
plotsol('h+l','pt0'); pause; plotsol('s+l','pt10'); pause; plotsol('sn1','pt10'); 