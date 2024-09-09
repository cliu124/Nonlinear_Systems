%% cmc spherical caps, basic script
close all; keep pphome; 
%% cmc spherical caps, init; pars will be overwritten in spcapinit
nx=10; al=0; h0=0; v0=0; a0=0; par=[h0; v0; al; a0]; % initial pars 
p=spcapinit(nx,par); plotsol(p); p=setfn(p,'cap1'); p.sol.ds=0.01;
p.nc.dsmax=0.2; p.nc.usrlam=[2 4];  p.plot.bpcmp=1; 
p.nc.ilam=[2 1]; p.nc.nq=1; p.fuha.qf=@qV; p.fuha.qfder=@qVjac; % cont V, free H 
p.sw.jac=0; p.sw.qjac=1; % using numerical Jacs for G, and approximate q_u  
p=cont(p,10); % go 
%% example of mesh-adaption; loading soln useful for testing parameters 
% such as ngen (number of adaption loops) and sig (frac of triangles to refine)
p=loadp('cap1','pt10','cap1r'); p.nc.ngen=1; p.nc.sig=0.2; p=oomeshada(p); 
p=cont(p,20); % continue refined solution 
%% branch plot of H over V, and some solution plots 
c=[2 1]; mclf(5); plotbra('cap1r','pt40',5,c); xlabel('V'); ylabel('H'); box on; 
plotsol('cap1','pt0');pause; plotsol('cap1','pt20');pause; plotsol('cap1r','pt40'); 