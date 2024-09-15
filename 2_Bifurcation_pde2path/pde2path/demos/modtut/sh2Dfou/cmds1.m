%% SH2D on small square, via dct 
close all; keep pphome; 
%% init 
p=[]; par=[-0.05 0 -1]; % lam,quad,cubic 
lx=2*pi; ly=lx; nx=40; ny=nx; dir='tt1'; 
p=shinit(p,lx,ly,nx,ny,par); p=setfn(p,dir); 
p.nc.lammax=2; p.nc.dsmax=0.1; p.sw.verb=2; p.nc.neig=8; p.ps=2; 
p=setbel(p,0,1e-6,5,@lssgmres); p.sw.eigssol=0; p.sw.verb=2;
p.sw.bifcheck=2; p.nc.mu1=1; tic; p=cont(p,10);toc
%% 1BP double, for small nx,ny, solving the ABE sometimes fails; then use gentau below
aux=[]; aux.soltol=1e-10; aux.m=2;  
p0=cswibra(dir,'bpt1',aux); % cswibra, then reset some parameters
p0.nc.dsmax=0.5; p0.nc.dsmin=0.1; p0.nc.lammax=1; p0.sol.ds=0.1; p0.nc.tol=1e-4; 
%% select and cont. 
p=seltau(p0,1,'b1'); p.sw.verb=0; p.nc.tol=1e-8; 
p=setbel(p,0,1e-6,5,@lssgmres); p.sw.eigssol=3; p.sw.verb=0;
tic; p=cont(p,12); toc
%% 2nd brnch 
p=seltau(p0,2,'b2'); p.sw.verb=0; p=cont(p,4); p.nc.tol=1e-8; p=cont(p,8); 
%% BD plot 
f=3; c=5; figure(f); clf; plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab',10); plotbra('b2',f,c,'cl','r','lab',8); ylabel('max(u)'); 
%% soln plot 
plotsol('b1'); pause; plotsol('b2'); 