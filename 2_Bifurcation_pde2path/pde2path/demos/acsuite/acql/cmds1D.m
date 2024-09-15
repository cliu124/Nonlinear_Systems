close all; keep pphome; % command templates for ql-AC; run cell-by-cell; 
%% 1D init 
p=[]; lx=1; nx=100; dim=1; par=[0.25; 0.5; 1; -0.3; 0.2];  % c0, lam, ga, del, epsi 
p=acinit(p,lx,nx,par,dim); p=setfn(p,'1D'); p=findbif(p,2); 
%% BP1, bifurcate in both direction, fold on the supercrit. branch yield a sort unstable part
p=swibra('1D','bpt1','1D1a',0.1); p=cont(p,30);  
p=swibra('1D','bpt1','1D1b',-0.05); p.sw.verb=2; p.nc.dsmax=0.05; p.sw.bifcheck=2; 
p=cont(p,30); 
%% BP2, just one direction cause symmetric
p=swibra('1D','bpt2','1D2',0.1); p.nc.dsmax=0.2; p=cont(p,30); 
%% BD plot
figure(3); clf; pcmp=0; plotbra('1D',3,pcmp,'cl','b','lsw',0); 
plotbra('1D1a',3,pcmp,'cl','k','lab',28); 
plotbra('1D1b',3,pcmp,'cl','k','lab',24); 
plotbra('1D2',3,pcmp,'cl','r','lab',[5 14 25]); 
xlabel('\lambda'); 
%% soln plot 
plotsol('1D1a','pt30'); pause; plotsol('1D1b','pt30'); pause; plotsol('1D2','pt10'); 
%% auxiliary: plot the nonlinear diff. coeff 
x=-0.1:0.1:1.8;y=0.25-0.3*x+0.2*x.^2; figure(11); plot(x,y); axis([-0.1 1.5 0.1 0.3]); 
legend('c(u)'); xlabel('u'); 
%% jaccheck; rel error is about 0.002, but cont works  
[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); spy(Gd>e1/2); 
