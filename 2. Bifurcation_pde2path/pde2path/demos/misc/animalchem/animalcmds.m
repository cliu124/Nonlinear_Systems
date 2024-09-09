%% 'Animal'-chemotaxis model
% u_t=0.25 Lap u- lam*div(u*grad v)+r u(1-u)
% v_t=Lap v+(u/(1+u)-v)
%%
close all; %clear all; 
% Trivial branch with small dlammax and lam=11.8 to find bif-points 
p=[]; p=animalinit(p); screenlayout(p); p.nc.smod=0; p=cont(p); 
%% First bifurcating branch 
q=swibra('p','bpt1','q',-0.1); 
q.sol.xi=0.2; q.nc.nsteps=5; q.nc.dlammax=1; q=cont(q); 
%% Other direction
q.sol.ds=-0.2; q.nc.smod=25;q.nc.nsteps=100; q=cont(q);
%% Continue trivial branch to left for plotting
p0=[];p0=animalinit(p0); screenlayout(p); p0.nc.smod=0; p0.sol.ds=-0.5; 
p0.nc.lammin=7.5; p0.nc.dlammax=1; p0=cont(p0);
%% Bifurcation-Diagram 
figure(3); clf; plotbraf('p','pt40',3,1,'lw',4,'cl','b');
%plotbraf('p0',3,1,'lw',4,'cl','b');
plotbraf('q',3,1,'lw',4,'lab',10);
xlabel('\lambda');ylabel('||u_1-1||_{L^1}/|\Omega|');
%% Plots
plotsolf('q','pt25',1,1,2); axis image; colormap copper; 
plottauf('p','bpt1',7,1,2); axis image; colormap copper; 
plottauf('p','bpt2',8,1,2); axis image; colormap copper;


