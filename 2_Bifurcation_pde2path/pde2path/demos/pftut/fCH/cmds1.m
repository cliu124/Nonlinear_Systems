close all; keep pphome; p=[];
%% fCH, init straight interface 
lx=2; yfac=1.5; ly=yfac*lx; p.pp=3; p.del=0; sw=1; wi=0.15; a=1.25; nx=45; 
ny=round(yfac*nx/2); p.fuha.wfu=@wfu2; p.qf=pi/2*lx; 
eta1=1; eta2=2; eps=0.25; mi=0; s=0; ga=0; sdir='i1'; 
par=[eta1; ga; eps; mi; eta2; s]; p=fchinits(p,lx,ly,nx,ny,par,wi,a,sw); 
p.nc.lammax=1;  p.nc.lammin=-1; p=setfn(p,sdir); p.sw.verb=1; p.fsol.fsol=0; 
stansavefu(p); p.sw.verb=2; p.plot.pstyle=1; p.nc.bisecmax=5; p.nc.neig=10; 
%[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); spy(Gd>e1/2); 
%% run mesh-ada,  OOPDE, yielding np\approx 5500 (see bottom for trulle) 
p.nc.ngen=1; p.nc.sig=0.75; p.sw.trul=0; p.fuha.e2rs=@e2rs_ad_hoc; p=oomeshada(p); p.np, stansavefu(p); 
%% now continue, first load and plot again (to check which mesh) 
p=loadp('i1','pt0','s1'); plotsol(p); p.np, pause 
p.plot.pstyle=2; p.nc.dsmax=0.1;  p.sol.ds=-1e-2;  
p.nc.ilam=[4 2 6]; % cont in [m,ga,s], m=mass, ga=mass-constrainnt, s=speed  
p.sw.spcalc=1; p.nc.amod=0;  p=cont(p,50); 
%% put BPs here for inspection 
p=swibra('s1','bpt4','du',0.01);
%% now the relevant branches 
p=swibra('s1','bpt3','p1',0.01); pause;  p=cont(p,20); 
p=swibra('p1','bpt1','p1-1',-0.05); pause; p=cont(p,20);  % secondary 
p=swibra('s1','bpt2','m1',0.01); pause;  p=cont(p,20); 
p=swibra('s1','bpt1','m2',0.01); pause;  p=cont(p,20); 
p=swibra('m1','bpt1','m1-1'); pause; p=cont(p,10);  % secondary, with drift 
%% plot BD, first at low m, gamma=3 
figure(3);clf;cmp=8; plotbra('s1','pt20',3,cmp,'cl','k','bplab',[1 2 3 4]); 
plotbra('p1','pt10',3,cmp,'cl','b','lab',10); 
plotbra('m1','pt40',3,cmp,'cl',p2pc('g1'),'lab',[40]); 
plotbra('m2','pt20',3,cmp,'cl',p2pc('g2'),'lab',20); 
xlabel('m'); ylabel('max(u)'); 
%% BD at larger m, cmp 8=max (good to see spots!) 
figure(3);clf;cmp=8; plotbra('s1','pt20',3,cmp,'cl','k','lab',0,'fp',1); 
plotbra('s1b','pt10',3,cmp,'cl','k','lab',0,'fp',5); 
plotbra('p1','pt10',3,cmp,'cl',p2pc('b1'),'lab',10); 
plotbra('p1-1','pt10',3,cmp,'cl','r','lab',5); xlabel('m'); ylabel('max(u)'); 
%% solution plots 
ps=2; myax='image'; v=[0,90]; 
plotsol('s1','pt0',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('p1','pt10',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('m1','pt40',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('m2','pt20',10,1,ps); axis(myax); view(v); nolti; pause 
%%
plotsol('s1','pt0',1,1,1); view([10,60]); axis([0 lx -ly -2.5]); nola 
%% tangent plot via calling swibra again, then some settings
q=swibra('s1','bpt2','du',0.01); nolti; view(0,90); axis image; 
title('\tau_1 at bpt2','FontSize',12); set(gca,'FontSize',14);
%% *********************************** experiment with trulle, 
% nice coarsening, but symmetry loss makes continuation less robust
p=loadp('s1','pt10','i2'); 
%% repeat for further reduction 
p.nc.ngen=2; p.nc.imax=20; op=troptions2D(); 
op.innerit=1; op.verbose=2; op.ppar=2; op.setids=@setidssq; 
op.Lup=2; op.Llow=0.2; p.nc.tol=1e-10; 
op.etafu=@etafua2D; p.trop=op;  % put options in p 
p.sw.trul=1; p.sw.ips=0; % interpolation switch, deals with boundaries 
p=oomeshada(p); p.np, plotsol(p,1,1,1); stansavefu(p);