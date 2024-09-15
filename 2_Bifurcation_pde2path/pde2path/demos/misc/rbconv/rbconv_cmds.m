%% Rayleigh-Benard convection demo
% Vorticity formulation with no-slip and stress-free boundary conditions on
% a rectangle.
%%
% Use the following as command templates and run cell-by-cell 
close all; keep pphome; 
%% Initialize with no-slip b.c. and find 2 bifurcation points on trivial branch 
ns=[]; ns=rbconv_init(ns); ns.sol.ds=10; ns.nc.nsteps=50; 
%% Set ds,nsteps outside of init since specific for findbif
ns=findbif(ns,2); 
%% Switch to continuation of bifurcating branches, use pmcont for speedup
ns1=swibra('ns','bpt1','ns1'); ns1.sol.ds=2; ns1=pmcont(ns1);
ns2=swibra('ns','bpt2','ns2'); ns2.sol.ds=2; ns2=pmcont(ns2);
%% Plot Bifurcation diagram 
figure(3);clf(3);cmp=0;
plotbraf('ns',3,cmp,'cl','b');
plotbraf('ns1',3,cmp,'cl','r','lab',9);
plotbraf('ns2',3,cmp,'cl','k','lab',7);
axis([739 940 0 3]); xlabel('R'); ylabel('max\psi');
%% Plot solutions 
p1=loadp('ns1','pt9'); arrowplot(p1,8,1); 
p2=loadp('ns2','pt7'); arrowplot(p2,9,1); 
%% Switch to stress-free boundary conditions and find bifurcations on trivial branch 
sf=[];sf=rbconv_init_stressfree(sf);
sf.sol.ds=10; sf.nc.nsteps=50; sf=findbif(sf,2);
%% First bifurcating branch (takes a while)
sf1=swibra('sf','bpt1','sf1',2); sf1.sol.dsmax=10; sf1.nc.nsteps=10; sf1=pmcont(sf1); 
%% Second bifurcating branch 
sf2=swibra('sf','bpt2','sf2',2); sf2.sol.dsmax=10; sf2=pmcont(sf2); 
%% Switch to (wrongly) detected 'bifurcating' branch
p=swibra('sf2','bpt1','sf2b'); p.nc.nsteps=10; p=cont(p); p=pmcont(p);
%% Switch in the opposite direction
% Fails: it recovers the first branch. 
p=swibra('sf2','bpt1','sf2b2');
p.sol.ds=-p.sol.ds; p.nc.nsteps=10; p=cont(p);
% This suggest a wrongly detected bifurcation point due to large stepsize.
%% The bifurcation is in fact 'imperfect'
% Slow backward continuation reveals this. However, using p.mesh.sympoi=1 in
% initialization the bifurcation is indeed a pitchfork.
p=loadp('sf2','pt10','sf3'); p=resetc(p);
p.sol.ds=-p.sol.ds; p.nc.dsmax=2; p.nc.nsteps=30; p=cont(p); p=pmcont(p);
%% Plot Bifurcation Diagram 
cmp=0;
plotbra('sf',3,cmp,'lw',5,'cl','b');
plotbra('sf1',3,cmp,'cl','r'); 
plotbra('sf2',3,cmp,'cl','k');
plotbra('sf2b',3,cmp,'cl','k');
plotbra('sf3',3,cmp,'cl','g');
axis([660 920 0 3]); xlabel('R'); ylabel('max\psi');
%% Plot solutions 
p1=loadp('sf1'); p2=loadp('sf2'); 
p3=loadp('sf2b'); p4=loadp('sf3'); 
arrowplot(p1,8,1); arrowplot(p2,9,1); 
arrowplot(p3,10,1); arrowplot(p4,11,1); 
