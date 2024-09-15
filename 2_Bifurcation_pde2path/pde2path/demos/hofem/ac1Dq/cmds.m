%% ac1D with quadr FEM 
close all; keep pphome; 
%% init and cont trivial branch 
par=[1 -0.1 1 0 0.5]; % c, lam, quintic, bdval, drift 
nx=20; p=acinit(2*pi,nx,par); p=setfn(p,'tr'); 
p.nc.dsmax=0.05; p.sw.bifcheck=2; mclf(10); spy(p.mat.K); p.nc.neig=4; p=cont(p,30);
%% switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); pause; p=cont(p,20);  
p=swibra('tr','bpt2','b2',0.1); pause; p=cont(p,20); 
p=swibra('tr','bpt3','b3',0.1); pause; p=cont(p,20); 
%%  BD plot
f=3; c=0; figure(f); clf;
plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','r','lab',10); 
plotbra('b2',f,c,'cl','b', 'lab',10); 
plotbra('b3',f,c,'cl','m','lsw',0); ylabel('||u||_2'); 
%% solution plots 
plotsol('b1','pt20'); nola; pause; plotsol('b2','pt20'); nola; pause; plotsol('b3','pt20'); nola
%% some mesh-refinement, first for one fixed solution
p=loadp('b1','pt20','b1');
p.nc.maxt=100; p.nc.dxmax=0.5; p.nc.redf=3; p.nc.dxmax=2; 
p=oomeshada(p,'ngen',1,'sig',0.8); plotsol(p,1,1,5);
mclf(10); spy(p.mat.K); 
%% continuation with mesh-adaption each amod-th step 
p=swibra('tr','bpt1','b1ref',0.1); p.nc.sig=0.5; p.nc.dxmax=1; 
p.nc.redf=p.hofem.femorder+1; p.nc.ngen=2; p.nc.sig=0.8; 
p.sw.errcheck=0; p.nc.amod=5; p.plot.pstyle=5; p.nc.lammax=2;  p=cont(p,20);  
figure(3); clf; plotbra(p,3,-1,'lsw',0,'labi',3); % plot error-est on branch