%% Allen-Cahn with global coupling; 
%
% G(u)=-c*lap(u)-lam*u-u^3+ga*u^5+del*<u^j>
%
% requires some global fields for Sherman Morrison (SM). See gclss and gcblss 
% Eigenvalues via mod of spcalc, which (essentially) uses eigs with gclss as solver 
% averages computed as <h(u)>=p.avvec*h(u); <h'(u) v>=avec*v, 
% where avec is computed (and passed on globally to gclss etc) in sGjac
close all; keep pphome; global p2pglob; 
%% First branch: brute=1: full Jac; SM faster for np=1000 (say) 
dim=2; lx=pi/2; nx=20; brute=0; 
par=[0.5; 0; 1; 0.5;  4]; % [d lam ga del j];   
p=[]; p=acgcinit(p,par,lx,nx,dim,brute); p.nc.dsmax=0.1; 
p=setfn(p,'2D0'); p=cont(p,40);  
%% bif branches 
p=swibra('2D0','bpt1','1',0.1); p.nc.dsmax=0.5; p=cont(p,10); 
p=swibra('2D0','bpt2','2',0.1); p.nc.dsmax=0.5; p=cont(p,10);
%% bra-plot 
figure(3); clf; pcmp=0;  plotbra('2D0','pt40',3,pcmp); 
plotbra('1','pt10',3,pcmp,'cl','b'); plotbra('2','pt10',3,pcmp,'cl','r'); 
%% soln plots
plotsol('1','pt10'); pause; plotsol('2','pt10'); 