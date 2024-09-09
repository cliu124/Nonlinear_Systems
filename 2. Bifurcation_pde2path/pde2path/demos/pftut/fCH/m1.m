close all; keep pphome; p=[];
%% fCH, commands for curved (initial) interfaces 
lx=1; yfac=2; ly=yfac*lx; p.pp=3; p.del=0; sw=1; wi=0.1; a=1.5; nx=50; 
lx=1; yfac=2; ly=yfac*lx; p.pp=3; p.del=0; sw=1; wi=0.125; a=1.5; nx=30; 
lx=1; yfac=3; ly=yfac*lx; p.pp=3; p.del=0; sw=1; wi=0.125; a=1.5; nx=30; 
ny=round(yfac*nx/2); p.fuha.wfu=@wfu2; 
eta1=1; eta2=2; eps=0.25; mi=-0*0.88; s=0; 
ga=0; sdir='c6i'; % nloop goes to pearling
par=[eta1; ga; eps; mi; eta2; s]; p=fchinits(p,lx,ly,nx,ny,par,wi,a,sw); 
p.nc.lammax=1;  p.nc.lammin=-1; p=setfn(p,sdir); p.sw.verb=1; p.fsol.fsol=0; 
stansavefu(p); p.sw.verb=2; p.plot.pstyle=1; p.nc.bisecmax=5; 
%[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); spy(Gd>e1/2); 
%% first Newton 
[u1,r1,i1,Gu,Glam,p]=nloop(p,p.u); fprintf('res=%g, iter=%g\n',r1,i1); 
p.u=u1; plotsol(p,1,1,p.plot.pstyle); stansavefu(p); p.sw.verb=2; 
%% run mesh-ada, first OOPDE, yielding np\approx 4160
p=loadp('c6i','pt0','c6a'); p.sw.trul=0; 
%%
p.nc.ngen=1; p.nc.sig=0.7; p=oomeshada(p); p.np 
stansavefu(p); 
%% trulle, yielding np\approx 5460
p=loadp('c6a','pt1'); p.np
%%
p.nc.ngen=2; p.nc.imax=20; op=troptions2D(); 
op.innerit=2; op.verbose=2; op.ppar=2; op.setids=@setidssq; 
p.nc.tol=1e-10; op.etafu=@etafua2D3; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; % interpolation switch, deals with boundaries 
p=oomeshada(p); p.np, plotsol(p,1,1,1); 
%%
p=setfn(p,'c6a'); stansavefu(p);
%% now continue 
p=loadp('c6a','pt1','s1'); plotsol(p); p.np
%%
p.nc.lammin=0; p.plot.pstyle=2; p.nc.tol=1e-8; p.nc.dsmax=0.5; 
p.nc.ilam=[4 2 6]; p.nc.lammin=-1; p.sol.ds=-1e-2;  
p.pm.mst=6; p.pm.resfac=1e-3; p.sw.bifcheck=1; p.nc.mu2=0.01; 
p.nc.foldtol=0.02; %p.sol.xi=1e-3; p.sol.xiq=1e-3;
p.sw.spcalc=1; p.nc.amod=0;  p=cont(p,50); 
%%
p=loadp('s1','pt5','s1b'); p.sol.ds=0.1; p.nc.dsmax=0.1; p.sw.bifcheck=2; p=cont(p,10); 

%%
p=swibra('s1','bpt8','du',0.1);
%%
p=swibra('s1','bpt1','q1',0.1); pause;  p=cont(p,15); 
%%
p=swibra('q1','bpt1','q1-1'); pause; p=cont(p,10);  
%%
p=swibra('s1','bpt9','q3',-0.1); pause; p=cont(p,10); 
%%
p=swibra('s1','bpt5','q5',0.1); p.plot.pstyle=2; pause; p.sw.bifcheck=0; p=cont(p,10); 
%%
p=swibra('q2','bpt1','q2-1'); pause; p=cont(p,10);  
%% plot BD, gamma=3, max=8 (good to see spots!) 
figure(3);clf;cmp=3; 
plotbra('s1','pt10',3,cmp,'cl','k','lab',0,'lp',7); 
plotbra('s1b','pt10',3,cmp,'cl','k','lab',0); 
plotbra('q1','pt10',3,cmp,'cl',p2pc('g1'),'lab',10); 
plotbra('q1-1','pt10',3,cmp,'cl',p2pc('g2'),'lab',10); 
xlabel('m'); ylabel('\gamma'); 
%%
figure(3);clf;cmp=3; 
plotbra('s1','pt20',3,cmp,'cl','k','lab',0,'fp',5); 
plotbra('q2','pt15',3,cmp,'cl',p2pc('b1'),'lab',10); 
plotbra('q2-1','pt10',3,cmp,'cl','r','lab',10); 
xlabel('m'); ylabel('\gamma'); 
%% solution plots 
ps=1; myax='tight'; v=[-40 60]; ps=2; myax='image'; v=[0,90]; %myax='tight'; 
plotsol('s1','bpt1',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('q1','pt10',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('q1-1','pt10',10,1,ps); axis(myax); view(v); nolti; 
%%
plotsol('q2','pt10',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('q2-1','pt10',10,1,ps); axis(myax); view(v); nolti; 
%% tangent plot via calling swibra again, then some settings
q=swibra('s1','bpt10','du',0.01); nolti; axis image; 
title('\tau_1 at bpt10','FontSize',p.plot.fs); set(gca,'FontSize',q.plot.fs);
%%
p=loadp('s1','bpt6','t1'); p.u(1:p.np)=p.u(1:p.np); 
t1=0; ts=[]; nc=0; dt=0.01; nt=3000; pmod=100; smod=100; p.mat.Kadv=0; 
%% the tint loop, repeat this cell until residual is small (here just once) 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalft); axis image; nolti


%% trulle, yielding np\approx 5460
p.nc.ngen=2; p.nc.imax=20; op=troptions2D(); 
op.innerit=1; op.verbose=2; op.ppar=2; op.setids=@setidssq; 
op.Lup=2; op.Llow=0.2; p.nc.tol=1e-10; 
op.etafu=@etafua2D; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; % interpolation switch, deals with boundaries 
p=oomeshada(p); p.np, plotsol(p,1,1,1); stansavefu(p);
%% now continue 
p=loadp('c2','pt0'); p.sol.ds=0.01; p.pm.resfac=1e-4; p=pmcont(p,150); 
%% other direction 
p=loadp('c2','pt5','c2b'); p.sol.ds=-0.1; p=pmcont(p,150);  
%% plot BD, gamma 
figure(3);clf;cmp=3; plotbra('c2','pt200',3,cmp,'cl','b','lab',[0 200]); 
plotbra('c2b','pt120',3,cmp,'cl',p2pc('b2'),'lab',[100]); 
xlabel('\eta_1'); ylabel('\gamma');title(['mass=-0.925']);
%% solution plots 
ps=1; myax='tight'; v=[-40 80]; ps=2; myax='image'; v=[0,90]; 
plotsol('c2','pt0',1,1,ps); nolti; axis(myax); view(v); pause 
ps=2; myax='image'; v=[0,90]; 
plotsol('c2','pt200',1,1,ps); nolti; axis(myax); view(v); pause 
plotsol('c2b','pt100',1,1,ps); nolti; axis(myax); view(v); 
%% check stab with tint
%p=loadp('c4','pt1','c4t'); 
t1=0; ts=[]; nc=0; dt=0.01; nt=2000; pmod=100; smod=100; p.mat.Kadv=0; 
%% the tint loop, repeat this cell until residual is small (here just once) 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalft); 
%% checking LSS 
p=loadp('c2','pt0','l2');  p.sw.verb=2; p.sw.spcalc=0; p.sw.para=2; p.sw.bifcheck=0; % 12s
p=setilup(p,1e-4); % 30s, 1 prec
%% 
p.nc.lammax=2.5; tic; p=cont(p,10); toc
%% cont in eps (only works till eps\approx 0.49 !!!) 
p=swiparf('c1tr','pt0','c1eps',[3 2]); p.nc.lammin=0; p.sol.ds=-0.05; 
p.plot.pstyle=1; p=pmcont(p,150);