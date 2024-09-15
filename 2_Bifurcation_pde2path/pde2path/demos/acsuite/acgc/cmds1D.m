%% Allen-Cahn with global coupling; 
%
% G(u)=-c*lap(u)-lam*u-u^3+ga*u^5+del*<u^j>
%
% requires some global fields for Sherman Morrison (SM). See gclss and gcblss 
% Eigenvalues via mod of spcalc, which (essentially) uses eigs with gclss as solver 
% averages computed as <h(u)>=p.avvec*h(u); <h'(u) v>=avec*v, 
% where avec is computed (and passed on globally to gclss etc) in sGjac
close all; keep pphome; global p2pglob; 
%% trivial branch: brute=1: full Jac; SM faster for np=1000 (say) 
dim=1; lx=pi/2; nx=30; icsw=0; brute=0; 
par=[0.5; 0; 1; 0.5;  4]; % [d lam ga del j];   
p=[]; p=acgcinit(p,par,lx,nx,dim,brute); p.nc.dsmax=0.1; 
huclean(p); p=setfn(p,'1D0'); p=cont(p,40);  
%% first two branches 
p=swibra('1D0','bpt1','1D1',0.1); p.nc.dsmax=0.5; p=cont(p,10); 
p=swibra('1D0','bpt2','1D2',0.1); p.nc.dsmax=0.5; p=cont(p,10);
%% bra-plot 
figure(3); clf; pcmp=0; plotbra('1D0','pt30',3,pcmp); 
plotbra('1D1','pt10',3,pcmp,'cl','b'); 
plotbra('1D2','pt10',3,pcmp,'cl','r'); 
%% soln plots
plotsol('1D1','pt10'); pause;  plotsol('1D2','pt10'); 
%% to switch on SM if it wasn't at startup
p=swibra('1D0','bpt1','1'); p.nc.dsmax=0.5; p.sw.bifcheck=0; 
p.jfac=0; p.fuha.lss=@gclss; p.fuha.blss=@gcblss; 
tic; p=cont(p,10); toc