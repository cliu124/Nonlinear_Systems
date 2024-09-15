%% GP1D, critical case sig=2
close all; keep pphome; 
%% init and cont of trivial branch, sig=2=critical  
p=[]; lam=1; s=1.8; sig=2; ga=-2; par=[lam s sig ga]; 
p=gpinit(p,10,400,par); p=setfn(p,'2/tr'); p.nc.neig=20; p.nc.lammax=5; 
p.sol.ds=-0.1; p=cont(p,10);
%% switch to first bifurcating branch and continue
p=swibra('2/tr','bpt1','2/b1',0.1);  p.file.smod=2; p.nc.amod=5; p=cont(p,40); 
%% 2ndary bif 
p=swibra('2/b1','bpt1','2/1a',0.1);  p.nc.amod=10; p.nc.sig=0.5; 
p.nc.lammax=10; p.fuha.outfu=@gpbra3; p.nc.dxmax=0.1;  p=cont(p,80); 
%% bifurcation diagram plotting 
f=3; c=6; ax=[0 5 0 2.2]; l='N'; %c=7; ax=[0 3.1 -15 1]; l='H'; 
figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('2/tr','pt10',f,c,'cl',[0.5 0.5 0.5],'lsw',0); 
plotbra('2/b1',f,c,'cl','k','lwst',2,'lab',38); 
plotbra('2/1a',f,c,'cl','b','lwst',2,'lab',26);
xlabel('\lambda'); ylabel(l); axis(ax); 
%% solution plots 
plotsol('2/b1','bpt1'); nola; grid on; pause;
plotsol('2/b1','pt38'); nola; grid on; pause; plotsol('2/1a','pt26'); nola; grid on
%% plot difference to pure NLS soliton 
f=3; c=7; clf(3); 
plotbra('2/1a','pt53',f,c,'cl','b','lwst',2,'lab',28,'fp',1);
xlabel('\lambda'); ylabel('||u(\cdot)-u_{NLS}(\cdot-s)||_2'); 
%% complex form for spectral problem.  
dir='2/b1'; p=real2c(dir,'pt16'); 
%% just spec for other solns 
dir='2/1a'; p=real2c(dir,'pt10'); 
%% -------- time integration, starting near b1/pt8, preparation -----------
p=loadp('1/b1','pt8','1-b1-8t'); p.prep=1; lam=p.u(p.nu+1); 
T=50*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.005; dtp=1; nc=0; 
N=1000;  nplot=floor(T/dtp); nt=ceil(T/dt)
%% the tint-loop  (can be repeated) 
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,0.01); 
%% plot tint data 
vv=[0,90]; fssplot('1-b1-8t',0,nplot,1,1,vv,'b1/pt8, \delta=0.01');
%% starting near 1a/pt9,
p=loadp('1/1a','pt9','1-1a-9t'); p.prep=1; lam=p.u(p.nu+1); 
T=90*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.005; dtp=1; nc=0; 
nplot=floor(T/dtp); N=1000; nt=ceil(T/dt)
%%
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,0.01); 
%%
vv=[0,90]; fssplot('1-1a-9t',0,nplot,1,1,vv,'1a/pt9, \delta=0.01');