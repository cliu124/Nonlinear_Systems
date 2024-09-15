%% ac1D with Chebychev-differentiation, NBC, 
close all; keep pphome; 
%% init 
p=[]; par=[1 -0.1 0 -1]; % c,lam,quad,cubic 
p=acinit(p,2,40,par); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=0.2; p.sw.verb=2; p.nc.neig=10; p.nc.bisecmax=10; 
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8; p.nc.eigref=-2; p=cont(p,30);
%% switch to first 4 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); pause; p=cont(p);  % spat.hom 
p=swibra('tr','bpt2','b2',0.1); pause; p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); pause; p=cont(p); 
p=swibra('tr','bpt4','b4',0.1); pause; p=cont(p); 
%% BD plot 
f=3; c=5; figure(f); clf; plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lsw',10);  plotbra('b2',f,c,'cl','r', 'lab',10); 
plotbra('b3',f,c,'cl','m', 'lab',5); plotbra('b4',f,c,'cl',p2pc('o1'), 'lab',5);
ylabel('max(u)'); axis([-0.1 7 0 1.7]); 
%% soln plots 
plotsol('b2','pt10'); pause; plotsol('b3','pt5'); pause; plotsol('b4','pt5'); 