%% ac1D, NBC, using dct 
close all; keep pphome; 
%% init 
p=[]; par=[0.5 -0.2 0 -1]; % c,lam,quad,cubic, 
p=acinit(p,2*pi,30,par); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=2; p.nc.dsmax=0.1; p.sw.verb=2; p.nc.neig=10;
p.sw.bifcheck=2; p.nc.mu1=1; p=cont(p,20);
%% switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p);  
p=swibra('tr','bpt2','b2',0.1); p=cont(p);
p=swibra('tr','bpt3','b3',0.1);  p=cont(p); 
%% 
f=3; c=5; figure(f); clf;
plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab',10); 
plotbra('b2',f,c,'cl','r', 'lab',10); 
plotbra('b3',f,c,'cl','m','lab',10); 
ylabel('max(u)'); 
%% soln plots 
plotsol('b2','pt10'); pause; plotsol('b3','pt10'); 
%% check diff-mat, full, but diag-dom
figure(10); clf; plotmat(p.mat.L); pause; clf; spy(p.mat.L); 