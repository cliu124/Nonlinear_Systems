%% SH1D, NBC, via dct, u in x-space 
close all; keep pphome; 
%% init 
p=[]; par=[-0.02 2 -1]; % lam,quad,cubic 
p=shinit(p,12*pi,100,par); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=2; p.nc.dsmax=0.01; p.sw.verb=2; p.nc.neig=10;
p.sw.bifcheck=2; p.nc.mu1=1; p=cont(p,8);
%% switch to first bifurcating branch and continue
p=swibra('tr','bpt1','b1',0.1); pause; p.nc.dsmax=0.02; p=cont(p,10); p.nc.dsmax=0.1; p=cont(p,30);  
%% snake
p=swibra('b1','bpt1','c1',0.1);pause; p.nc.dsmax=0.1;  p=cont(p,95);
%% BD plot 
f=3; c=3; figure(f); clf;
plotbra('tr','pt10',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab',32); 
plotbra('c1',f,c,'cl','r', 'lab',[40 80],'lp',102); 
ylabel('||u||_2'); 
%% soln plots 
plotsol('b1','pt30'); pause; plotsol('c1','pt40'); pause; plotsol('c1','pt80'); 
%% inspect diff.matrix 
figure(10); clf; plotmat(p.mat.L); 