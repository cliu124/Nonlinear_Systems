close all; keep pphome; 
%% commands for Schnakenberg on a small 2D domain with hex lattice
p=[]; kc=sqrt(sqrt(2)-1); lx=2*pi/kc; ly=lx/sqrt(3); par=[3.3, 0, 60]; 
nx=30; ny=round(nx/sqrt(3)); p=schnakinit(p,lx,ly,nx,ny,par); % init 
p.pm.resfac=1e-4; p.sol.ds=-0.1; p=setfn(p,'hom');plotsol(p); p.sw.verb=2; p.sw.bifcheck=2; pause 
p.sw.verb=1; p.np, p=findbif(p,2); p=cont(p,20); 
%% hex via qswibra, continue in both directions
p0=qswibra('hom','bpt1'); p0.nc.dsmin=0.1;p0.nc.bisecmax=4; 
p0=setbelilup(p0,1,1e-4,5,1e-6,200); p0.nc.mu2=0.02;  
%% hex
p=seltau(p0,3,'h-',2); p.sol.ds=0.1; p=pmcont(p,15);
p=seltau(p0,3,'h+',2); p.sol.ds=-0.1; p=pmcont(p,20);
%% s+ and s- via gentau
p=gentau(p0,1,'s+'); p.sol.ds=-0.05; p=cont(p,20); 
p=gentau(p0,1,'s-'); p.file.smod=5; p=cont(p,20); 
%% b+ and b-
p=swibra('s+','bpt1','b+',-0.05); p.nc.dsmin=0.05; p.pm.resfac=1e-3; p.file.smod=5; p=pmcont(p,15);
%p=swibra('s-','bpt1','b-',0.05); p.nc.dsmin=0.05; p.pm.resfac=1e-3; p.file.smod=5; p=pmcont(p,10);
%% plot BD 
fnr=3; figure(fnr); clf; plotbra('hom',fnr,4,'cl','k','lsw',0);
plotbra('s+',fnr,4,'cl','b','lab',15); plotbra('s-',fnr,5,'cl','b','lsw',0); 
plotbra('h+',fnr,4, 'cl','m','lab',10);  plotbra('h-',fnr,5, 'cl','m','lab',10); 
plotbra('b+',fnr,4,'cl','r','lab',5,'lp',19); 
plotbra('b-',fnr,5,'cl',p2pc('o1'),'lab',5); 
axis([2.4 3.25 0.7 5.8]); xlabel('\lambda'); ylabel('max/min(u_1)'); 
%% soln plot  
plotsol('h+','pt10');pause; plotsol('h-','pt10');pause; 
plotsol('b+','pt5'); pause; plotsol('b-','pt5'); pause; plotsol('s+','pt15'); 