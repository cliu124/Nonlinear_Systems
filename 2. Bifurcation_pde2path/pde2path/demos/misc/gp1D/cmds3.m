%% GP, supercritical case sig=3
close all; keep pphome; 
%% init and cont of trivial branch, sig=3, supercrit
p=[]; lam=1; s=1; del=1; sig=3; ga=-2; par=[lam s sig ga]; 
p=gpinit(p,10,1000,par); p=setfn(p,'3/tr'); p.nc.neig=20; p.nc.mu1=1; p.nc.mu2=0.1; 
p.sol.ds=-0.1; p.nc.maxt=2000; p=cont(p,10); 
%% switch to first 2 bifurcating branches and continue
p=swibra('3/tr','bpt1','3/b1',0.01); p.nc.amod=5; p.nc.dsmax=0.04; p.file.smod=1; p=cont(p,35);
p=swibra('3/tr','bpt2','3/b2',0.1); p.nc.amod=5; p=cont(p,30); 
%% 2ndary bif 
p=swibra('3/b1','bpt1','3/1a',-0.1); pause; p=pmcont(p,50); 
%% bifurcation diagram plotting 
f=11; c=6; ax=[0.5 2.3 0 2]; l='N'; %c=7; ax=[0 3.1 -15 1]; l='H'; 
figure(f); clf; plotbra('3/tr','pt10',f,c,'cl',[0.5 0.5 0.5],'lsw',0); 
plotbra('3/b1','pt50',f,c,'cl','k','lab',42); 
plotbra('3/1a','pt50',f,c,'cl','b','lwst',2); 
xlabel('\lambda'); ylabel(l); axis(ax); 
%% solution plots 
plotsol('3/b1','bpt1'); nola; pause; plotsol('3/b1','pt42'); nola; pause; 
plotsol('3/1a','pt34');
%% complex form for spectral problem.  
dir='3/b1'; p=real2c(dir,'pt4',1,1e-4); %p=real2c(dir,'bpt1',2,1e-4); 
%% do cont backward to check stab at low ampl. 
clf(2); p=resetc(p); p.sol.ds=0.01; p.nc.dsmax=0.05; p=cont(p,5);
%% complex form for spectral problem.  
dir='3/b1'; p=real2c(dir,'pt42'); 
%% complex form for spectral problem.  
dir='3/1a'; p=real2c(dir,'pt48'); pause; p=resetc(p); p.sol.ds=0.05; p.nc.dsmax=0.05; p=cont(p,20);
%% more careful cont near max of b1 branch 
p=real2c('3/b1','pt10'); p.sol.dsmax=0.05; p.sw.verb=2; p=setfn(p,'3/b1b'); 
p.fuha.outfu=@gpbra2; p=cont(p,5);  
%% -------- time integration, starting near b1/bpt1, preparation 
p=loadp('3/b1','bpt1','3-b1-bp1'); p.prep=1; lam=p.u(p.nu+1); 
T=3*pi/lam; tp=2*pi/lam; dtp=1; 
t1=0; ts=[]; dt=0.001; dtp=0.05; nc=0; N=1000; nt=ceil(T/dt)
%% excess mass (delta>0) yields blowup 
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,0.001); 
figure(10); clf; plot(ts(1,:),ts(2,:),ts(1,:),ts(3,:)); 
%%
vv=[120,40]; fssplot('3-b1-bp1',0,75,1,1,vv,'bpt1, \delta=0.001');
%% stable for mass defect 
p=loadp('3/b1','bpt1','3-b1-bp1b'); p.prep=1; lam=p.u(p.nu+1); 
T=10*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.005; dtp=0.5; 
nc=0; N=1000; nplot=floor(T/dtp); nt=ceil(T/dt); nf=20; 
%%
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,-0.001,nf); 
%%
figure(10); clf; plot(ts(1,:),ts(2,:),ts(1,:),ts(3,:)); 
vv=[0,90]; fssplot('3-b1-bp1b',0,2*nplot,1,1,vv,'bpt1, \delta=-0.001');
%% tint from point from unstable branch with lower N yiels trans.to 'small lam' GS 
p=loadp('3/b1','pt34','3-b1-34t'); p.prep=1; lam=p.u(p.nu+1); 
T=10*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.01; dtp=0.5; 
nc=0; N=1000; nplot=floor(T/dtp); nt=ceil(T/dt)
%%
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,-0.01,nf); 
%%
figure(10); clf; plot(ts(1,:),ts(2,:)); %,ts(1,:),ts(3,:)); 
vv=[0,90]; fssplot('3-b1-34t',0,nplot,1,1,vv,'b1/pt42, \delta=-0.01');
%% load point from unstable branch with lower N, conv.to 'small lam' GS 
p=loadp('3/1a','pt34', '3-1a-34t'); p.prep=1; lam=p.u(p.nu+1); 
T=10*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.01; dtp=0.5; 
nc=0; N=1000; nplot=floor(T/dtp); nt=ceil(T/dt), nf=20; 
%%
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,0.00,nf); 
%%
figure(10); clf; plot(ts(1,:),ts(2,:),ts(1,:),ts(3,:));    
vv=[0,90]; clf(2); fssplot('3-1a-34t',0,nplot,1,1,vv,'1a/pt48, \delta=0');

