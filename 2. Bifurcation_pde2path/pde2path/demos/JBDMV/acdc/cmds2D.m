%% 2D dead core 
close all; keep pphome; 
%% init, 
p=[]; par=[-0.1 0.85 0.05 1 0.3 8 0.6]; dim=2; sw=1; del=1.2; 
p=acinit(p,1,120,par,dim,sw,del); p=setfn(p,'2'); p.Jdel=1e-8; plotsol(p); 
p.nc.ngen=1;  p.nc.sig=0.6; p.nc.maxt=12000; % mesh-adapt.settings 
p.nc.lammax=40; p.nc.tol=1e-6; % slightly relax tol 
%% first 4 BPs, then cont with larger ds 
p=findbif(p,4);  p.nc.dsmax=5; p.nc.dlammax=2; p=cont(p,20); 
%% center DC branch, with meshada
p=swibra('2','bpt1','c1',-0.1); p.nc.amod=10; p=cont(p,50);  % center DC
%% 2nd branch 
p=swibra('2','bpt2','c2',0.01); p.nc.sig=0.8; p.nc.amod=10; p=cont(p,100); 
%% BD, L^2
f=3; c=0; figure(f); clf;
plotbra('2',f,c,'cl','k'); plotbra('c1',f,c,'cl','b','lab',[49]); 
plotbra('c2',f,c,'cl','m','lab',[100]); ylabel('||u||_2'); 
%% sol plots; 
v=[0 90]; 
mypsol2D('c1','pt49',v); pause; mypsol2D('c2','pt100',v); pause; 
plotsol('c2','pt100',1,1,3);  nola; colormap parula; 