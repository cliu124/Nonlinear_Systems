%% extended Brusselator, 1D, basics 
close all; keep pphome; 
%% init, and cont of hom.steady branch 
ndim=1; dir='hom1d'; p=[]; lx=pi/1.4; nx=100; 
par=[0.95; 2.75; 1; 1; 0.01; 0.1; 1]; % a, b, c, d, Du, Dv, Dw 
p=bruinit(p,lx,nx,par,ndim); p=setfn(p,dir); p.sw.bifcheck=2;
p=initeig(p,4); p.nc.neig=[5, 5]; % init omv (compute guesses for eval shifts) 
p.sw.verb=2; p.nc.mu2=5e-3; p=cont(p,20); 
%% 2 steady bifs from BPs to Turing branches 
p=swibra('hom1d','bpt1','1ds1',0.01); p.nc.dsmax=0.05; p=cont(p,20);
p=swibra('hom1d','bpt2','1ds2',0.01); p.nc.dsmax=0.05; p=cont(p,20);
%% Hopf bifs from HBPs
para=4; ds=0.1; dsmax=0.2; xi=1e-3; nsteps=15; figure(2); clf; 
aux=[]; aux.tl=30; 
for bnr=1:3
switch bnr
 case 1; p=hoswibra('hom1d','hpt1',ds,para,'1dh1',aux); 
 case 2; p=hoswibra('hom1d','hpt2',ds,para,'1dh2',aux); p.hopf.smod=1; 
 case 3; p=hoswibra('hom1d','hpt3',ds,para,'1dh3',aux); 
end 
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.sw.verb=0; 
p=setbel(p,2,1e-4,5,@lss); t1=tic; p=cont(p,nsteps); toc(t1) 
end 
%% 2ndary Bifs  
para=4; ds=0.2; p.nc.dsmax=0.2; aux=[]; aux.tl=20; p=hoswibra('1ds1','hpt1',ds,para,'1ds1h1',aux); 
p.hopf.xi=1e-2; p.nc.dsmax=1; p.sw.verb=1; p=setbel(p,1,1e-4,5,@lss);
tic; p=cont(p,15); toc
% p=hoswibra('1ds2','hpt1',ds,para,'1ds2h1',aux); % quite similar 
%% plot BD, L^2
bpcmp=9; wnr=3; figure(wnr); clf
plotbra('hom1d','bpt3',3,bpcmp,'cl','k'); 
plotbra('1ds1','pt20',3,bpcmp,'cl','b'); 
plotbra('1ds2',3,bpcmp,'cl',[0 0.25 1]); 
plotbra('1dh1',3,bpcmp,'cl',[1 0.5 0],'lab',10); 
plotbra('1dh2',3,bpcmp,'cl',[1 0.5 0],'lab',[5,10]); 
plotbra('1dh3',3,bpcmp,'cl',[1 0.5 0],'lab',[5,10]); 
plotbra('1ds1h1','pt15',3,bpcmp,'lab',[10],'cl','r'); 
xlabel('b'); ylabel('max(u)');
%% plot solns
plotsolf('1ds1','hpt1',1,1,1); figure(1); 
title('s1/hpt1'); xlabel('x'); ylabel('u_1'); pause; 
hoplotf('1dh1','pt10',1,1); figure(1); title('h1/pt10'); colormap cool; 
xlabel('x'); ylabel('t'); zlabel('u_1');  v=[30 60]; view(v); pause
hoplotf('1dh3','pt5',1,1); figure(1); title('h3/pt5'); colormap cool;  
xlabel('x'); ylabel('t'); zlabel('u_1');  view(v); pause
hoplotf('1ds1h1','pt10',1,1); figure(1); title('s1h1/pt10'); colormap cool; 
xlabel('x'); ylabel('t'); zlabel('u_1');  view(v);
%% stab check, init 
p=loadp('1ds1h1','pt10'); hoplot(p,1,1); dir='stab1d'; 
p.u(1:p.nu)=p.hopf.y(1:p.nu,1); u0=p.u(1:p.nu); p=setfn(p,dir);
ts=[]; t0=0; npp=50; pmod=50; smod=5; tsmod=1; nc=0; 
%% tint
nt=1500; [p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); 
figure(4); clf; plot(ts(1,:), ts(2,:)); % point-vals
figure(5); clf; plot(ts(1,:), ts(3,:)); % difference in norm
set(gca,'FontSize',p.plot.fs); axis tight; 
xlabel('t'); ylabel('||u(t)-u_0||_{\infty}'); 
%% x-t plot: (see in ts if there's something interesting after np periods, then 
%               plot around there)
dir='stab1d'; wnr=2; cmp=1; incr=5; timesl=3; 
switch timesl % choose timelslice to plot
    case 0; si=0; nt=1*npp/smod; zti=[0.8 1.2]; yti=[0 5];  
    case 1; si=350; nt=2*npp/smod; zti=[1 1.5]; yti=[5 10 15 20 25 30 35 40 45 50]; 
    case 2; si=1300; nt=2*npp/smod; zti=[0.6 1.2 1.8]; yti=[150 155 160 165 170]; 
end
vv=[30,70]; tintplot1d(dir,si,incr,nt,wnr,cmp,vv);
colormap cool; shading interp; view(30,70); grid off;
set(gca,'YTick',yti); set(gca,'ZTick',zti);
xlabel('x'); ylabel('t'); zlabel('u_1')
%% (compute) and plot multipliers
aux.nfloq=100; muv1=floqap('1dh3','pt10',aux); 
%% check hploc, here near eigref(2)
p=loadp('hom1d','pt1','hploc_test'); p=hploc(p,2); 
para=4; ds=0.1; dsmax=0.2; xi=1e-3; figure(2); clf; aux.tl=30; 
p=hoswibra('hploc_test','hpt1',ds,para,'1dh1b',aux); 
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.sw.verb=0; p=cont(p,5);