%% demo fhneps, clear,  - init and findbif 
close all; format compact; keep pphome; p=[]; p=FHNinit(p,10,50); 
%pars: 1=eps; 2=vel; others of coupling function g=p_3+p_4*u+p_5*u^2+p_6*u^3
par=[6; 0; 0; 0; 0; 0;]; u=zeros(p.np,1); v=u; p.u=[u;v; par]; p=findbif(p,1); 
%% cell 1 - swibra to stationary front for one step
p=swibra('init','bpt1','prep/swibra',0.1); p.nc.lammin=0.5; p=cont(p,1);
%% cell 2 - add velocity and phase equation
p=swiparf('prep/swibra','pt1','prep/decoup',[1;2]); 
p.nc.nq=1; p.nc.xiq=0.1; p.fuha.qf=@qf; p.fuha.qfder=@qfder;
p.sol.ds=-0.1; p.nc.lammin=0.1; p.usrlam=[0.4,0.3,0.2]; p.sw.bifcheck=0; 
p.nc.dsmax=0.5; p.nc.dlammax=0.5; p.sw.bprint=2; clf(2); p=cont(p,30); 
%% cell 3 - refine mesh
p=loadp('prep/decoup','pt20','prep/decoup'); 
p=meshada(p,'ngen',10,'sig',1e-4);  p.fuha.savefu(p); 
%% cell 4 - turn on cubic coefficient
p=swiparf('prep/decoup','pt21','prep/cubic',[6;2]); 
p.nc.lammin=-1; p.nc.lammax=1; 
p.usrlam=[0.25,0.5,0.75]; p.sol.ds=-0.1; clf(2); p=cont(p,20); 
%% cell 5 - find bifurcation point
p=swiparf('prep/cubic','pt5','stand',[4;2]); p.sw.bifcheck=1; 
p.nc.lammin=-2; p.nc.lammax=1.5; p.sol.ds=0.1; 
p.usrlam=[]; clf(2); p=cont(p,25); plotsol(p,6);
%% cell 6 - switch to branch with nonzero velocity
p=swibra('stand','bpt1','travel',0.1); p.nc.dsmax=0.1; 
p.plot.bpcmp=2; clf(2); p=cont(p,20); 
%%
p=swibra('stand','bpt1','travel-',-0.01); p.nc.dsmax=0.1; 
p.plot.bpcmp=2; clf(2); p=cont(p,20); 