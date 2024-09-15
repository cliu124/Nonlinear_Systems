%% example commands for mesh-modification, 2D SH
keep pphome; close all; p=[]; 
%% init and zero-branch with a very coarse mesh 
p=[]; lx=2*pi; nx=round(3*lx); ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=1.3; par=[lam; nu];  
sw.ref=0; sw.sym=1; % ref=#mesh-refs, sym=1: cc-mesh 
p=shinit(p,nx,lx,ly,ndim,par,sw); huclean(p); p=setfn(p,'hc'); p.plot.pstyle=1; 
p.np, p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.05; 
p.file.smod=1; p.sw.bifcheck=2; p=cont(p,10); 
%% swibra to hex/stripes 
aux=[]; aux.m=2;   p0=qswibra('hc','bpt1',aux); p0.sw.bifcheck=1; 
%% hexagons 
p=seltau(p0,2,'hc1',2); p.sol.ds=-0.05; p=pmcont(p,40); 
%% stripes 
p=gentau(p0,[0 1],'sc1'); p=pmcont(p,20); 
%% soln plots 
plotsol('hc1','pt10',1,1,2);
%% settings for trullekrul
p=loadp('hc1','pt10','h1r'); p=resetc(p); p.plot.shsw=0; 
p.plot.pstyle=1; plotsol(p); p.np, 
op=troptions2D(); % load default trullerup-options, then overload some 
op.innerit=2; op.verbose=2; op.ppar=2; op.setids=@setidssq;  op.Lup=2; op.Llow=0.075;
op.etafu=@etafua2D; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; p.nc.ngen=1; p=oomeshada(p); plotsol(p);  stansavefu(p);
%% meshada each 3rd step; here p.sw.ips=2 seems inportant, otherwise extrapol. errors! 
p=loadp('h1r','pt0','h1ada'); p.nc.amod=3; p.nc.ngen=1; p.sw.ips=2; p=cont(p,10); 
%% use pa2pb to move soln to new mesh, then cont 
q=[]; lx=2*pi; nx=round(5*lx); ly=lx/sqrt(3); ndim=2; lam=-0.001; nu=1.3; par=[lam; nu];  
sw.ref=0; sw.sym=1; % ref=#mesh-refs, sym=1: cc-mesh 
q=shinit(p,nx,lx,ly,ndim,par,sw); huclean(q); q=setfn(q,'dummy'); p.plot.pstyle=2; 
p=loadp('hc1','pt20','h1r'); 
%%
q=pa2pb(p,q); 
%%
q=cont(q,10); 
