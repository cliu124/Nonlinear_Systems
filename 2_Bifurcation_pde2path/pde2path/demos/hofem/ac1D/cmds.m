%% ac1d with Lobatto FEM 
close all; keep pphome; 
%% init, with femorder=m, and cont 
par=[1 -0.1 1 0 0]; nx=6; m=4; p=acinit(2*pi,nx,par,m); p=setfn(p,'tr'); 
p.nc.dsmax=0.05; p.sw.bifcheck=2; mclf(10); spy(p.mat.K); p.nc.neig=4; 
bpoints=[1/16, 1/4, 9/16]; b=0*bpoints;  p=cont(p);
%% switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); b(1)=getlam(p); p=cont(p,20);  
p=swibra('tr','bpt2','b2',0.1); b(2)=getlam(p); p=cont(p,20); 
p=swibra('tr','bpt3','b3',0.1); b(3)=getlam(p); p=cont(p,20); 
%% error in BPs: P1: 0.0103951, P2: 0.000349952, m=3: 0.000565, m=4: 0.000215185 
for i=1:3; fprintf('%g ',bpoints(i)-b(i)); end; % error in BPs 
fprintf(',   ||bpoints-b||_2=%g \n', norm(bpoints-b,2));   
%%  BD plot
f=3; c=0; figure(f); clf;
plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','r','lab',10); 
plotbra('b2',f,c,'cl','b', 'lab',10); 
plotbra('b3',f,c,'cl','m','lsw',0); ylabel('||u||_2'); 
%% solution plots 
plotsol('b1','pt10'); nola; yticks([0.4 0.8]); pause; plotsol('b2','pt10');
%% some mesh-refinement, first for one fixed solution
p=loadp('b1','pt20','b1'); % p.hofem.femorder=3; % can also change order here! 
p.nc.maxt=100; p.nc.dxmax=0.5; p.nc.redf=p.hofem.femorder; p.nc.dxmax=2; 
p=oomeshada(p,'ngen',1,'sig',0.5); plotsol(p,1,1,5);
mclf(10); spy(p.mat.K); 
%% continuation with mesh-adaption each amod-th step 
p=swibra('tr','bpt1','b1ref',0.1); p.nc.sig=0.5; p.nc.dxmax=1; p.hofem.femorder=3; 
p.nc.redf=p.hofem.femorder; p.nc.ngen=2; p.nc.sig=0.3; 
p.sw.errcheck=0; p.nc.amod=5; p.plot.pstyle=5; p.nc.lammax=2;  p=cont(p,20);  
figure(3); clf; plotbra(p,3,-1,'lsw',0,'labi',3); % plot error-est on branch