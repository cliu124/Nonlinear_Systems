%% AC on graph 
close all; keep pphome; 
%% init and cont trivial branch 
p=[]; par=[1 -0.2 0 -1]; % parameters c,lam,quad,cubic 
sw=4; np=250; sw=0; np='G2'; % uncomment 2nd part to load fixed graph 
p=acinit(p,par,np,sw); p=setfn(p,'tr2'); p.nc.dsmax=0.1;  
p.nc.neig=50; p.nc.nsteps=1000; p=findbif(p,5); % run continuation 
%% swibra 
p=swibra('tr2','bpt1','c1',0.1); p=cont(p,10);  
p=swibra('tr2','bpt2','c2',0.1); pause; p=cont(p,10);  
p=swibra('tr2','bpt3','c3',0.1); pause; p=cont(p,10);  
p=swibra('tr2','bpt3','c5',0.1); pause; p=cont(p,10);  
%% BD plot
f=3; c=5; figure(f); clf;
plotbra('tr2',f,c,'cl','k','lsw',0); plotbra('c1',f,c,'cl','b','lsw',0); 
plotbra('c2',f,c,'cl','r', 'lab',10); plotbra('c3',f,c,'cl','m','lab',10); 
ylabel('max(u)'); 
%% solution plot
plotsol('c2','pt10'); pause; plotsol('c3','pt10'); 