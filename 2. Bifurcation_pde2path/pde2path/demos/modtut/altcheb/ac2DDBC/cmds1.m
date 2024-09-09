%% ac2D, DBC, chebychev-setup 
close all; keep pphome; 
%% init 
p=[]; par=[0.5 -0.1 0 -1]; % c,lam,quad,cubic
lx=2; ly=1; nx=40; ny=20; p=acinit(p,lx,ly,nx,ny,par); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=1; p.sw.verb=2; p.nc.neig=10;
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8; p=cont(p,15);
%% switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); pause; p=cont(p);  
p=swibra('tr','bpt2','b2',0.1); pause; p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); p=cont(p); 
%% BD plot 
f=3; c=5; figure(f); clf;
plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab',10); 
plotbra('b2',f,c,'cl','r', 'lab',10); 
plotbra('b3',f,c,'cl','m','lab',10); 
ylabel('max(u)'); axis([0 11 0 3.5]); 
%% soln plots 
plotsol('b1','pt10'); pause; plotsol('b2','pt10'); pause; plotsol('b3','pt10'); 
%%
figure(10); clf; spy(p.mat.L); 