%% ac2D, NBC, chebychev-setup 
close all; keep pphome; 
%% init 
p=[]; par=[0.5 -0.1 0 -1]; % c,lam,quad,cubic
lx=2; ly=1; nx=40; ny=20; Lfn='Ls'; % filename for diff.op.matrix L 
p=acinit(p,lx,ly,nx,ny,par,Lfn); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=1; p.sw.verb=2; p.nc.neig=10; p.sw.jac=1; 
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8;  p.sw.verb=2; 
%p=findbif(p,4);
%% switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); p=cont(p,10);  
p=swibra('tr','bpt2','b2',0.1); pause; p=cont(p,10); 
p=swibra('tr','bpt3','b3',0.1); pause; p=cont(p,10); 
p=swibra('tr','bpt4','b4',0.1); pause; p=cont(p,10); 
%% BD plot 
f=3; c=5; figure(f); clf;
plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lsw',0); 
plotbra('b2',f,c,'cl','r', 'lab',10); 
plotbra('b3',f,c,'cl','m','lab',5); 
plotbra('b4',f,c,'cl',p2pc('o1'),'lab',5); 
ylabel('max(u)'); axis([0 4.5 0 2.2]);
%% soln plots 
plotsol('b2','pt10'); pause; plotsol('b3','pt5'); pause; plotsol('b4','pt5'); 
%%
spy(p.mat.L); 
%%
[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); em=max(max(Gd)); figure(11); spy(Gd>em/10); 