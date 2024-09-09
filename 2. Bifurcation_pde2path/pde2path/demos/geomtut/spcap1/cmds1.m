%% cmc spherical cap, cont in (H,V)
close all; keep pphome;  global p2pglob; 
p2pglob.tsw=1; p2pglob.vi=[20,40]; p2pglob.edc='k'; % plotting controls
%% init; pars will be overwritten in scinit
nx=12; h0=0; v0=0; a0=0; par=[h0; v0; a0]; % initial pars 
p=scinit(nx,par); p=setfn(p,'cap1'); p.sol.ds=0.1; 
p.sw.jac=0;  % numerical (0) or functional (1) jacs for G, speed no problem 
p.sw.qjac=1; % numerical (0, too slow), hence functional (1) derivative for q
p.nc.ilam=[2 1]; p.nc.nq=1; p.fuha.qf=@qfV; p.fuha.qfder=@qjacV; % cont in V, 
p.plot.bpcmp=1; p.nc.usrlam=[2 4]; % cmp for branch-plot, vals for forced output 
p=cont(p,5); % go 
%% alternate cont and mesh refinement based on triangle areas 
p=loadp('cap1','pt5','cap1r'); p.sw.nobdref=0; p.sw.rlong=1; p.file.smod=2; 
sig=0.2; for i=1:10; p=refineX(p,sig); p=cont(p,5); end 
%% just cont in H (no constraints) 
p=loadp('cap1','pt0','Hcont');p.nc.nq=0;p.nc.ilam=1;p.sol.ds=-0.1;p.plot.bpcmp=5; 
sig=0.2; for i=1:4; p=cont(p,5); p=refineX(p,sig); end % alternate ref. and cont 
%% branch plot of H over V and A over V, 
c=[5 1]; mclf(5); plotbra('cap1r','pt30',5,c); xlabel('V'); ylabel('H'); box on; 
c=[6 1]; mclf(7); plotbra('cap1r','pt30',7,c); xlabel('A'); ylabel('H'); box on; 
%% some solution plots 
p2pglob.tsw=3; plotsol('cap1','pt0'); pause; plotsol('cap1','pt5'); pause; 
plotsol('cap1r','pt18'); 
%% HK plots 
p2pglob.edc='k'; p2pglob.cb=1; plotHK('cap1r','pt40'); 
%% cont in H, A;  same branch as for H,V  
nx=12; h0=0; v0=0; a0=0; par=[h0; v0; a0]; % initial pars 
p=scinit(nx,par); p=setfn(p,'cap1A'); plotsol(p); p.sol.ds=-0.1;p.nc.dsmax=0.2; 
p.nc.ilam=[1 3]; p.nc.nq=1; p.fuha.qf=@qfA; p.sw.qjac=1; p.fuha.qfder=@qjacA;  
huclean(p); p.plot.bpcmp=3; p2pglob.tsw=2; % use (A,H) as title 
p.file.smod=2; p.sw.jac=0; p=cont(p,10); % go 
%% mesh-ada, based on triangle areas 
p=loadp('cap1A','pt10','cap1Ar'); sig=0.1; p.sw.rlong=1; p=refineX(p,sig); p=cont(p,2);  
%% cont further; with repeated mesh adapation; 
% or repeat p=cont(p,2); p=refineX(p,0.05); from the command line 
sig=0.1; for i=1:8;  p=cont(p,2); p=refineX(p,sig); end  
p=cont(p,4); % some further cont
%% branch plot of A over H, 
c=[1 5]; mclf(5); plotbra('cap1Ar','pt30',5,c,'lab',[10, 30]); 
box on; xlabel('H'); ylabel('A');
%% some solution plots 
plotsol('cap1A','pt10'); pause; plotsol('cap1Ar','pt30');
%% branch plot of V over A, and some solution plots 
c=[6 5]; mclf(5); plotbra('cap1Ar','pt28',5,c); xlabel('A'); ylabel('V');