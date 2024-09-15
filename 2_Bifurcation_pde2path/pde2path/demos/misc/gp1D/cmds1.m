%% GP1D, subcritical case sig=1
close all; keep pphome; 
%% init and cont of trivial branch, sig=1=subcritical 
p=[]; lam=1; s=2; sig=1; ga=-2; par=[lam s sig ga]; 
p=gpinit(p,10,300,par); p=setfn(p,'tr'); p.nc.neig=20; 
p.sol.ds=-0.1; p=cont(p,10);
%% switch to first 2 bifurcating branches and continue
p=swibra('tr','bpt1','1/b1',-0.01);  p.file.smod=2; p.nc.amod=5; p=cont(p,30); 
p=swibra('tr','bpt2','1/b2',0.1);  p.nc.amod=5; p=cont(p,30); 
%% 2ndary bif 
p=swibra('1/b1','bpt1','1/1a',0.1); p=cont(p,30); 
%% deflation; if it fails, try replacing cosh(x+s) below by cosh(x-s) 
p=loadp('1/1a','pt14'); plotsol(p); nu=p.nu; pause
x=getpte(p); x=x'; global p2pglob; al3=1; p=deflinit(p,al3); % init deflation solver 
p.defl.nsw=2; % reset selected deflation parameters/switches, here norm 
p.nc.tol=1e-4; p.defl.al3=1; p.sw.newt=0; % 0,1=nloop, 2=NLEQ1, 3=fsolve
ug=p.u; amp=max(p.u(1:nu)); s=p.u(p.nu+2); 
ug(1:nu)=ug(1:nu)+amp./cosh(2*(x+s)); plotsolu(p,ug,1,1,1);
p.nc.tol=1e-5; [u, p]=deflsol(p,ug); plotsol(p); 
p=postdefl(p,2,'1/d1a',-0.1);  p.nc.sig=0.8; 
p.nc.tol=1e-8; p.sw.newt=1; p.nc.amod=10; p.nc.dxmax=0.1; p=cont(p,30);
%%
p=loadp('1/d1a','pt4','1/d1b'); p.sol.ds=-p.sol.ds; p=cont(p,20); % other direction 
%% bifurcation diagram plotting 
f=3; c=6; ax=[0 3.1 0 7]; l='N'; %c=7; ax=[0 3.1 -15 1]; l='H'; 
figure(f); clf; % f=figure-Nr, c=component number (of branch) 
plotbra('tr','pt10',f,c,'cl',[0.5 0.5 0.5],'lsw',0); 
plotbra('1/b1',f,c,'cl','k','lab',[4,8]); 
plotbra('1/b2',f,c,'cl','m','lab',11); 
plotbra('1/1a',f,c,'cl','b','lab',9);
plotbra('1/d1a',f,c,'cl',p2pc('o1'),'lab',21); 
plotbra('1/d1b',f,c,'cl',p2pc('o2'),'lab',18); 
xlabel('\lambda'); ylabel(l); axis(ax); 
%% solution plots 
plotsol('1/b1','pt4'); nola; pause; plotsol('1/b1','pt8'); nola; pause;
plotsol('1/1a','pt9'); nola; pause; plotsol('1/b2','pt11'); nola; pause; 
plotsol('1/d1a','pt21'); nola; pause; plotsol('1/d1b','pt18'); nola; 
%% complex form for spectral problem.  
dir='1/b1'; p=real2c(dir,'pt8'); 
%% do cont backward to check stab at low ampl. 
clf(2); p=resetc(p); p.sol.ds=-0.01; p.nc.dsmax=0.05; p=cont(p,13);
%%
plotspec('1/b1pt8v','pt2',6);pause; plotspec('1/b1pt8v','pt10',6);
%% just spec for other solns 
dir='1/b2'; p=real2c(dir,'pt11'); pause; dir='1/1a'; p=real2c(dir,'pt9'); 
%% branches from deflation, may need NLEQ1 (newt=2) (and small tol) for init 
newt=1; dir='1/d1a'; p=real2c(dir,'pt21',newt,1e-2); pause; dir='1/d1b'; 
p=real2c(dir,'pt18',newt,1e-2); 
%% check if cont also works 
clf(2); p.sol.ds=0.1; p.nc.dsmax=0.1; p=cont(p,10); 
%% -----------   BP-cont in s (separation)  ------------------------
p=spcontini('1/b1','bpt1',2,'1/bpc1'); % init BP cont with s new prim. par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; p.sol.ds=-0.05;        
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
p.fuha.outfu=@stanbra; p.nc.amod=0; p.sw.bifcheck=0; p.nc.dsmax=0.2; 
%[Ja, Jn]=spjaccheck(p); pause % check impl. of spjac (comment out when fine) 
p=cont(p,50); 
%% plot BD and solns for BP-cont
clf(3); plotbra('1/bpc1','pt40',3,1,'lab',[30 40]); 
ylabel('\lambda'); %axis([0.65 2 0 8]); 
plotsol('1/bpc1','pt30'); nola; pause; plotsol('1/bpc1','pt40'); nola; 
%% -------- time integration, starting near b1/pt8, preparation -----------
p=loadp('1/b1','pt8','1-b1-8t'); p.prep=1; lam=p.u(p.nu+1); 
T=50*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.005; dtp=1; nc=0; 
N=1000;  nplot=floor(T/dtp); nt=ceil(T/dt), nf=0; 
%% the tint-loop  (can be repeated) 
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,0.01); 
%% plot tint data 
vv=[0,90]; fssplot('1-b1-8t',0,nplot,1,1,vv,'b1/pt8, \delta=0.01');
%% starting near 1a/pt9,
p=loadp('1/1a','pt9','1-1a-9t'); p.prep=1; lam=p.u(p.nu+1); 
T=90*pi/lam; tp=2*pi/lam; t1=0; ts=[]; dt=0.005; dtp=1; nc=0; 
nplot=floor(T/dtp); N=1000; nt=ceil(T/dt)
%%
[p,t1,ts,nc]=fsstint(p,t1,ts,dt,dtp,nt,nc,N,0.01,0); 
%%
vv=[0,90]; fssplot('1-1a-9t',0,nplot,1,1,vv,'1a/pt9, \delta=0.01');