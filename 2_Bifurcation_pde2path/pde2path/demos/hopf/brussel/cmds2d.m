%% extended Brusselator 2D, needs initeig to prepare HP detection 
close all; keep pphome; 
%% C1: init hom branch, with INITEIG, then use cont to find bifurcations 
Du=0.01; Dv=0.1; Dw=1; c=1; d=1; a=0.95; b=2.75; lx=pi/2;
ndim=2; dir='hom2d'; p=[]; nx=60; par=[a b c d Du Dv Dw]; 
p=bruinit(p,lx,nx,par,ndim); p=setfn(p,dir); p.sw.spcalc=0; p.nc.mu2=0.5e-2; 
p=initeig(p,4); p.nc.neig=[3, 3]; % init omv (compute guesses for eval shifts) 
p.sw.bifcheck=2;  p=cont(p,30); % cont with just 3 evals near 0 and near om1
%% C2: Bif from HPs, h1 
para=4; ds=0.1; dsmax=0.5; aux=[]; aux.tl=15; 
p=hoswibra('hom2d','hpt1',ds,para,'2dh1',aux); 
p.hopf.ax='unif'; % use same axis for all snapshots from Hopf orbit
p.hopf.fltol=1e-3; % poor accuracy of mu_1 due to coarse mesh in t! 
p.hopf.xi=1e-2; p.nc.dsmax=dsmax; 
AMG=1; p.sw.verb=3; p.hopf.flcheck=0; % set AMG=1 if ilupack available 
if ~AMG; p=setbel(p,1,1e-3,5,@lss); % bel significantly faster here 
else p=setilup(p,1e-3,50); p.fuha.blss=@lssAMG; end 
tic; p=cont(p,10); toc 
%% C3: steady bif to patterns; cont yields 2ndary Hopf
p=swibra('hom2d','bpt1','2ds1',0.01); p.sw.spcalc=1; p.nc.dsmax=0.01;
p.nc.dsmax=0.1; p=cont(p,2); % 2 initial steps 
p=meshada(p,'ngen',2); % then some mesh-adaption
p=pmcont(p,40); % use pmcont to avoid branch jumping
%% C4: 2ndary Hopf bifs from Turing branch 
aux=[]; aux.tl=15; ds=0.1; aux.tl=20; 
p=hoswibra('2ds1','hpt1',ds,para,'2ds1h1',aux); pause
p.file.smod=1; p.nc.dsmax=0.1; 
AMG=1; p.sw.verb=3; p.hopf.flcheck=0; 
if ~AMG; p=setbel(p,1,1e-3,5,@lss); % bel significantly faster here 
else p=setilup(p,1e-3,50); p.fuha.blss=@lssAMG; end 
p.nc.tol=1e-6; tic; p=cont(p,5); toc
%% C5: plot BD, L^2
bpcmp=9; wnr=3; figure(wnr); clf
plotbra('hom2d','hpt1',3,bpcmp,'cl','k'); 
plotbra('hom2d',3,bpcmp,'cl','k', 'lwst',2); 
plotbra('2dh1','pt10',3,bpcmp,'cl','r', 'lab',5); 
plotbra('2ds1','pt30',3,bpcmp,'cl','b'); 
plotbra('2ds1h1','pt10',3,bpcmp,'cl','m', 'lab',5,'lwst',2);
xlabel('b'); ylabel('max(u)');
axis([2.75 2.95 2.9 3.4]); set(gca,'XTick',[2.8 2.9]);
text(2.79,3.2,'s1h1','fontsize',16); text(2.84,3.2,'h1','fontsize',16);
text(2.76,3,'s1','fontsize',16);
%% C6: plot solns 
plotsolf('2ds1','hpt1',1,1,2); axis image; set(gca,'YTick',[-0.2 0.2]); 
aux.ytics=[]; aux.xtics=[]; aux.ztics=[1 1.5]; hoplotf('2ds1h1','pt5',4,1,aux);
%% C7: tint preparation 
p=loadp('2ds1h1','pt5'); hoplot(p,1,1); dir='stab1'; 
p.u(1:p.nu)=p.hopf.y(1:p.nu,1); u0=p.u(1:p.nu); p=setfn(p,dir);
ts=[]; t0=0; npp=50; nt=200; pmod=50; smod=5; tsmod=1; nc=0; 
%% C8: tint, repeat this cell for longer time integration 
nt=400; [p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); 
%% C9: plot time-series and solns 
figure(4); clf; plot(ts(1,:), ts(2,:)); % point-vals
figure(5); clf; plot(ts(1,:), ts(3,:)); % difference in norm
set(gca,'FontSize',p.plot.fs); axis tight; 
xlabel('t'); ylabel('||u(t)-u_0||_{\infty}'); 
dir='stab1'; si=0; incr=100; nt=1*npp/smod; wnr=2; cmp=1; pstyle=2; 
ind=si:incr:7*incr; lay=[4 2]; notic=1; 
tintplot2d(dir,ind,wnr,cmp,pstyle,lay,notic);
%% C10: movie 
p=loadp('2ds1h1','pt5');  mov=homov(p,1,1); mymov2avi(mov,'m1');
% avconv -i m1.avi -q 10 -vf 'setpts=8*PTS' m1.ogv 
%% C11: plot multipliers (somewhat slow!) 
aux.nfloq=400; muv1=floqap('2ds1h1','pt5',aux); 