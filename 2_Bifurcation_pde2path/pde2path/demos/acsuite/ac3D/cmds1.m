%% 3DAC, demo to check basics, on a rather coarse mesh
close all; keep pphome; 
%% C1: init, nx=30->np=10350
p=[]; par=[1 0.3 1 0]; lx=2*pi; ly=3*pi/2; lz=pi; nx=30; 
p=acinit(p,lx,ly,lz,nx,par); p.np, plotsol(p,1,1,1); p.fuha.spjac=@spjac; 
p.nc.ilam=2; p.nc.lammax=2; p.sol.ds=0.1; p.nc.dsmax=0.2; p=setfn(p,'tr'); p0=p; 
%% C2: use findbif to locate the BPs; here we comment out this cell and  
% alternatively use bplocin the next cell (more efficient for large nx) 
tic; p=findbif(p,3); toc
%% C3: localize BPs by setting lam to near a (here known) BP and using bploc 
%p.branch=[bradat(p); p.fuha.outfu(p,p.u)]; 
p=p0; p=cont(p,1); % 1 step, only to generate tangent 
tic; p=setlam(p,0.4); p=bploc(p); % localize 1st BP 
p=setlam(p,0.6); p=bploc(p); p=setlam(p,0.8); p=bploc(p); toc % 2nd and 3rd BP 
%% C4: switch to first 2 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.05); p=cont(p);
p=swibra('tr','bpt2','b2',0.05); p=cont(p); 
%% C5: bifurcation diagram plotting, 
f=3; c=0; figure(f); clf;
plotbra('tr', f,c,'cl','k'); plotbra('b1','lab',10,'cl','b'); 
plotbra('b2',f,c,'cl','r', 'lab',20); axis([0.2 2 0 28]); 
%% C6: solution plots 
plotsol('b1','pt10',1,1,1); pause; plotsol('b2','pt10',1,1,2,'alpha',0.5); axis image; 
%% C7: continue in some other param, here the BC coefficent d 
p=swiparf('b1','pt10','b1-dc',4); p.sol.ds=0.1; p.nc.lammin=-1; p.nc.lammax=2; 
p=resetc(p); clf(2); p.sw.eigssol=0;  p=cont(p,10); 
clf(1); ps=3; plotsol('b1-dc','pt10',1,1,ps); % solution plot 
%% C7b: 'cut-away'-plot
figure(1); clf; p=loadp('b1-dc','pt10'); 
p.plot.cut=[0 -10 -10]; % cutting at x<x1, y<y1, z<z1, here y1=0, x1,z1 outside dom. 
ps=4; plotsol(p,1,1,ps); % solution plot 
%% C8: fold-continuation 
p=spcontini('b1','fpt1',3,'b1f');   % init fold continuation with par 3 new primary par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.01; p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
tic; p=cont(p,10); toc