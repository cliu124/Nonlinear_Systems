%% commands for SH as consistent 2nd order sys, 1D  
keep pphome; close all; p=[]; 
%% init and zero-branch on long domain 
lx=20*pi; nx=round(30*lx); ly=0.1; ny=1; ndim=1; lam=-0.05; nu=2; 
par=[lam; nu]; p=shinit(p,nx,lx,ly,ndim,par); p.plot.pmod=10; 
p=setfn(p,'1D0b'); p=findbif(p,4);
%% 4 Turing-branches 
p=swibra('1D0b','bpt1','1D1b',0.01); p.sw.bifcheck=0; p=cont(p,40);
p=swibra('1D0b','bpt2','1D2b',0.01); p.sw.bifcheck=0; p=cont(p,40);
p=swibra('1D0b','bpt3','1D3b',0.01); p.sw.bifcheck=0; p=cont(p,40);
p=swibra('1D0b','bpt4','1D4b',0.01); p.sw.bifcheck=0; p=cont(p,40);
%% BD plot, H  
pcmp=6; fnr=4; figure(fnr); clf; 
plotbra('1D0b','bpt4',fnr,pcmp,'cl','k'); plotbra('1D1b','pt40',fnr,pcmp,'cl','k');
plotbra('1D2b','pt40',fnr,pcmp,'cl','b'); plotbra('1D3b','pt40',fnr,pcmp,'cl','r');
plotbra('1D4b','pt40',fnr,pcmp,'cl','m'); plot([-0.6 0.5],[0 0],'--k');  
axis([-0.55 0.1 -0.05 0.075]);  xlabel('\lambda'); ylabel('H'); 