%% ac1Dnlbc demo 
close all; keep pphome; 
%% init and cont in d=right BC 
c=1; lam=-0.1; ga=1; d=0; al=1; beta=0; 
p=[]; par=[c lam ga d al beta]; 
p=acinit(p,5,100,par); p=setfn(p,'a');
p.nc.ilam=4; p.nc.lammax=2; p.sol.ds=0.1; p.nc.dsmax=0.5; 
p.usrlam=[0.2 1 2]; p.file.smod=10; p.nc.dsmax=0.05; p.nc.neig=50; p.sw.verb=2; 
p.sw.foldcheck=1; p=cont(p,70);
%% branchplot and soln plot
figure(3); clf; plotbra('a','pt70',3,0); 
plotsol('a','pt6',1,1,1); pause; plotsol('a','pt35',1,1,1); pause; 
plotsol('a','pt64',1,1,1);
%% fold-continuation 
p=spcontini('a','fpt1',3,'af1');  % init fold continuation with par 3 new primary par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.01;                   % new stepsize in new primary parameter
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
p.file.smod=5; p=cont(p,10); 
%%
clf(3); plotbra('af1','pt20',3,4); % plot BD for fold-cont
%% switch back to regular continuation from one of the fold points
p=spcontexit('af1','pt10','a1-a'); 
p.nc.dsmax=0.5; p.sw.bifcheck=0; p.plot.bpcmp=0; p.sw.spcalc=1; 
p.sol.ds=-1e-3; clf(2); p=cont(p,1); p=cont(p,10); % continue in one direction 
%% cont in beta 
p=swiparf('a','pt64','b',6); p.nc.dsmax=0.1; p.usrlam=[0.2 1 2]; p=cont(p,30);
%%
figure(3); clf; plotbra('b',3,0);
%%
plotsol('b','pt4',1,1,1); pause; plotsol('b','pt13',1,1,1); pause; plotsol('b','pt24',1,1,1);