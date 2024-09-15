close all; keep pphome; p=[];
%% fCH, half circle 
lx=1.5*pi; yfac=1/2; ly=yfac*lx; p.pp=3; p.del=0; sw=2; wi=0.15; a=1.25; nx=130; 
ny=round(yfac*nx/2); p.fuha.wfu=@wfu2; p.qf=0.25*pi/lx; 
eta1=1; eta2=2; eps=0.25; mi=0; s=0; ga=0; sdir='i2'; 
par=[eta1; ga; eps; mi; eta2; s]; p=fchinits(p,lx,ly,nx,ny,par,wi,a,sw); 
p.nc.lammax=1;  p.nc.lammin=-1; p=setfn(p,sdir); p.sw.verb=1; p.fsol.fsol=0; 
stansavefu(p); p.sw.verb=2; p.plot.pstyle=1; p.nc.bisecmax=5; 
%[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); spy(Gd>e1/2); 
%% first Newton 
[u1,r1,i1,Gu,Glam,p]=nloop(p,p.u); fprintf('res=%g, iter=%g\n',r1,i1); 
p.u=u1; plotsol(p,1,1,p.plot.pstyle); stansavefu(p); p.sw.verb=2; 
%% run mesh-ada,  OOPDE, yielding np\approx 5500 (see bottom for trulle) 
p=loadp('i2','pt0','i2'); p.sw.trul=0; 
%% sig=0.5, then 0.2
p.nc.ngen=1; p.nc.sig=0.5; p=oomeshada(p); plotsol(p,1,1,1); p.np, stansavefu(p); 
%% now continue, first load and plot again, cause many different versions tried
%p=loadp('i2','pt0','c1'); plotsol(p); p.np, pause 
p=loadp('c1','pt15','c1'); plotsol(p); p.np, pause 
p.sw.bifcheck=2; 
p.plot.pstyle=2; p.nc.dsmax=0.01;  p.sol.ds=1e-2;  p.nc.neig=20;% p.sw.bifcheck=2; 
p.nc.ilam=[1 2 6]; % cont in m=mass, ga=mass-constrainnt, s=phase constraint 
p.nc.lammax=2; 
p.sw.spcalc=1; p.nc.amod=0;  p=pmcont(p,50); 
%% other direction 
p=loadp('s1','pt5','s1b'); p.sol.ds=0.01; p.nc.dsmax=0.1; p.sw.bifcheck=1; p=cont(p,10); 
%% put BPs here for inspection 
p=swibra('c1','bpt1','du',0.01);
%% now the relevant branches 
p=swibra('c1','bpt3','cp3',0.05); pause; %huclean(p); 
p.nc.neig=20; p.sw.bifcheck=2; p.sw.verb=2; p.nc.dsmax=0.05; 
p.nc.mu1=0.5; p.nc.mu2=0.025; 
p=cont(p,40); 
%p=swibra('c1','bpt2','cp2',0.05); pause; p.nc.neig=20; p=pmcont(p,20); 
%p=swibra('c1','bpt3','cp3',0.05); pause; p.nc.neig=20; p=pmcont(p,20); 
%p=swibra('p1','bpt1','p1-1'); pause; p=cont(p,20);  % secondary 
%%
p=swibra('s1b','bpt1','m1',0.01); pause;  p=cont(p,20); 
p=swibra('m1','bpt1','m1-1'); pause; p=cont(p,10);  % secondary, with drift 
%% plot BD, first at low m, gamma=3 
figure(3);clf;cmp=3; plotbra('c1','pt20',3,cmp,'cl','k'); 
plotbra('cp1','pt10',3,cmp,'cl','b','lab',0); 
plotbra('cp2','pt40',3,cmp,'cl','r','lab',0); 
plotbra('cp3','pt40',3,cmp,'cl','m','lab',0); 
%% solution plots 
ps=1; myax='tight'; v=[-40 60]; ps=2; myax='image'; v=[0,90]; 
plotsol('c1','bpt1',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('cp1','pt10',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('cp2','pt40',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('cp3','pt40',10,1,ps); axis(myax); view(v); nolti; 
%%
plotsol('m1','pt10',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('m1','pt20',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('m1-1','pt10',10,1,ps); axis(myax); view(v); nolti; 
%% tangent plot via calling swibra again, then some settings
q=swibra('s1b','bpt1','du',0.01); nolti; axis image; 
title('\tau_1 at bpt1','FontSize',14); set(gca,'FontSize',14);
%% time-integration 
p=loadp('s1','bpt8','t1'); p.u(1:p.np)=p.u(1:p.np); 
t1=0; ts=[]; nc=0; dt=0.01; nt=3000; pmod=100; smod=100; p.mat.Kadv=0; 
%% the tint loop, repeat this cell until residual is small (here just once) 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalft); axis image; nolti
%% run Newton on tint-soln to go to 'steady state' 
[u1,r1,i1,Gu,Glam,p]=nloop(p,p.u); fprintf('res=%g, iter=%g\n',r1,i1); 
p.u=u1; plotsol(p,1,1,p.plot.pstyle); 
p=setfn(p,'sp1a'); p=resetc(p); stansavefu(p); p.sw.verb=1; 
%% continue
p=loadp('sp1a','pt0'); p.sw.bifcheck=0; p.sw.spcalc=0; 
p=resetc(p); p.sol.restart=1; p.sol.ds=-0.01; p=cont(p,80);
%% other direction 
p=loadp('t1','pt0','t1'); p.sol.ds=-0.01; p.sw.bifcheck=0; p.sw.spcalc=0; 
p=cont(p,80); 
%% BD plot 
figure(3);clf;cmp=3;
plotbra('sp1a','p80',3,cmp,'cl','r','lab',[0 80]); 
plotbra('sp1b','p80',3,cmp,'cl','m','lab',80); xlabel('m'); ylabel('\gamma'); 
%% soln plot 
v=[90,90]; plotsol('sp1a','pt80',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('sp1b','pt80',10,1,ps); axis(myax); view(v); nolti; pause 
%% trulle, mainly for coarsening, but symmetry loss makes BDs less clear 
p=loadp('i2','pt1','i3'); 
%% repeat for further reduction 
p.nc.ngen=2; p.nc.imax=20; op=troptions2D(); 
op.innerit=1; op.verbose=2; op.ppar=2; op.setids=@setidssq; 
op.Lup=2; op.Llow=0.2; p.nc.tol=1e-10; 
op.etafu=@etafua2D; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; % interpolation switch, deals with boundaries 
p=oomeshada(p); p.np, plotsol(p,1,1,1); stansavefu(p);