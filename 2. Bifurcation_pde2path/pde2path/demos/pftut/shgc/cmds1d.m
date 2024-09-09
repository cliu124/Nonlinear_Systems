%% commands for SH with global coupling: globals avec, cvec and LU in p2pglob
% for gclss and lsslueigs will be put into p2pglob
keep pphome; close all; p=[]; global p2pglob; 
%% init and zero-branch 
lx=10*pi; nx=round(50*lx); ly=1; ndim=1; lam=-0.05; ga=2; nu=2; par=[lam; nu; ga]; 
p=shinit(p,nx,lx,ly,ndim,par); p=setfn(p,'1D0'); p=findbif(p,1);
%% first Turing-branch
p=swibra('1D0','bpt1','1D1',0.1); p.nc.dsmax=0.05; tic; p=cont(p,50); toc
%% localized slanted snake  
p=swibra('1D1','bpt2','1Ds1',0.01); p.nc.dsmax=0.1; p.sw.bifcheck=0; p=cont(p,110);
%% BD plot L2
pcmp=4; figure(3); clf; plotbra('1D1','pt44',3,pcmp,'cl','k');
plotbra('1Ds1','pt110',3,pcmp,'cl','b','lab',[10 50]);
axis([-0.5 1 0 0.8]); xlabel('\lambda'); ylabel('||u||'); 
%% soln plots 
plotsol('1Ds1','pt10'); xlabel(''); pause; plotsol('1Ds1','pt50'); xlabel('');