%% ac2Dpbc, here periodic in y 
close all; keep pphome; 
%% init, parameters [c lambda gamma d] 
p=[]; par=[1 -0.5 1 0]; lx=2*pi; ly=pi; p=acinit(p,lx,ly,40,par); 
p.nc.lammax=5; p.sol.ds=0.1; p.nc.dsmax=0.1; p=setfn(p,'tr'); 
%% first continue in d
p.nc.ilam=4; p=cont(p,4); 
%% find branchpoints for cont in lam  
p=resetc(p); p.nc.ilam=2; p=findbif(p,3); p=cont(p,30);
%% switch to first 3 bifurcating branches and continue; 
clf(2); 
for i=1:3; is=mat2str(i);     
    p=swibra('tr',['bpt' is],['b' is] ,0.05); p=cont(p,30); 
end
%% plot BD 
f=3; c=0; figure(f); clf;
plotbra('tr','pt50',f,c,'cl','k','lab',44,'bplab',1); 
plotbra('b1','pt30',f,c,'cl','b','lab',30); 
plotbra('b2','pt30',f,c,'cl','m','lsw',0); plotbra('b3','pt30',f,c,'cl','r','lsw',0); 
%% solution plots 
plotsol('tr','bpt1',1,1,3); pause; plotsol('tr','pt40',1,1,3); pause
plotsol('b1','pt30',1,1,3); 
%% fold-continuation 
p=spcontini('b1','fpt1',3,'b1f');   % init fold continuation with par 3 new primary par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.01; p.nc.tol=1e-5; p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
tic; p.nc.amod=0; p=cont(p,10); toc
%% switch back to regular continuation from one of the fold points
p=spcontexit('b1f','pt10','b1-a'); clf(2); p.sw.spcalc=1; 
p.nc.dsmax=0.5; p.sw.bifcheck=0;  p.nc.amod=0; p.plot.bpcmp=0; 
p.sol.ds=-1e-3; p=cont(p,1); p=cont(p,15); % continue in one direction 
p=loadp('b1-a','pt1','b1-b'); p.sol.ds=-p.sol.ds; % and the other one 
p.plot.bpcmp=0; p=cont(p,8); 
%% a cell to test mesh-adaption for pBC (2D) 
p=loadp('b2','pt10'); plotsol(p,1,1,1); view(-10,70); pause 
p.nc.ngen=1; p.nc.sig=0.3; p.nc.bddisty=0.4; %p.fuha.e2rs=@e2rs_ad_hoc; 
p=oomeshada(p); view(-10,60); pause 
p=cont(p,2); p=oomeshadac(p); p.nc.amod=10; p=cont(p,2); 