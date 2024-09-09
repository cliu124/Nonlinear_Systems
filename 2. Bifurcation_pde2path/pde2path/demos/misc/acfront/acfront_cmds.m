%% Allen-Cahn front demo
% 
%%
% Use the following as command templates and run cell-by-cell 
close all; keep pphome; 
%% Find bifurcation points from trivial branch 
p=[]; p=acfront_init(p); screenlayout(p); p=findbif(p,4); 
%% 
p=swibra('p','bpt1','q',0.2); p.sw.bifchecksw=0; p.sw.errchecksw=1; p.nc.amod=10; p=cont(p);
%% Add phase condition to simulate comoving frame for front-type solution
p=loadp('q','pt10','r');
p=resetc(p); p.nc.ilam=[2;3]; p.restart=1; % continue in symmetry breaking parameter and velocity
p.nc.nq=1; 
p.nc.amod=5; % mesh adaption all 5 steps
p.fuha.qf=@acfront_qf;   p.fuha.qfder=@acfront_qfder;
p.sw.qjac=1;
p.nc.lammax=20;
p.plot.bpcmp=3; clf(2); % bifurcation plots with parameter 3 as vertical axis
p=cont(p,20); % run continuation (there is a known explicit solution, 
%               velocity is linear in parameter 2)
