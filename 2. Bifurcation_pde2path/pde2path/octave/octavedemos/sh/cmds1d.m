%% commands for SH as consistent 2nd order sys, 1D  
keep pphome; close all; p=[]; 
%% init and zero-branch 
lx=6*pi; nx=round(40*lx); ly=0.1; ny=1; ndim=1; lam=-0.05; 
nu=2; % just some value for nu>nu0=sqrt(27/38)\approx 0.843
par=[lam; nu]; p=shinit(p,nx,lx,ly,ndim,par); huclean(p); 
p=setfn(p,'1D0'); p=findbif(p,4);
%% first Turing-branch
p=swibra('1D0','bpt1','1D1',0.01); p=cont(p,40);
%% two snaking branches (front and localized  
p=swibra('1D1','bpt1','1Ds1',0.01); p.nc.dsmax=0.1; p=cont(p,50);
p=swibra('1D1','bpt2','1Ds2',0.01); p.nc.dsmax=0.05; p=cont(p,70); 
%% BD plot L2
pcmp=3; figure(3); clf; plotbra('1D1','pt30',3,pcmp,'cl','k','lab',[30]);
plotbra('1Ds1','pt70',3,pcmp,'cl','b','lab',[80 90]);
plotbra('1Ds2','pt60',3,pcmp,'cl','r','lab',[80 100]);
xlabel('\lambda'); ylabel('||u||_2/sqrt(|\Omega|)'); 
%% plots (incl. Hamiltonian H and Fourier) 
hfplot('1D1','pt30'); pause; hfplot('1Ds1','pt50'); pause; 
hfplot('1Ds2','pt50'); pause; hfplot('1Ds2','pt30'); 
%% fold-cont 
p=spcontini('1D1','fpt1',2,'1D1f');   % init fold cont with par 2 new prim. par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new param. for plotting
p.sol.ds=0.1; p.nc.lammax=4; p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
[Ja, Jn]=spjaccheck(p); pause % check impl. of spjac (comment out when fine) 
tic; p=cont(p,10); toc
clf(3); plotbra('1D1f','pt10',3,1); % plot BD for fold-cont
%% switch back to regular continuation from one of the fold points
p=spcontexit('1D1f','pt10','1D1-a'); p.nc.dsmax=0.2; p.sw.bifcheck=0; 
p.plot.bpcmp=0; p.nc.lammin=-5; p.sol.ds=1e-2; clf(2); 
p=cont(p,1); p=cont(p,20); % continue in one direction, 1 initial step for saving
%% other direction 
p=loadp('1D1-a','pt1','1D1-b'); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p=cont(p,50); 
%% plot new branches
figure(3); clf; f=3; c=0; 
plotbra('1D1-a',f,c,'lsw',0); plotbra('1D1-b','pt50',f,c,'lsw',0); 
plotbra('1D1','pt40',f,c,'cl','r','lsw',0); 