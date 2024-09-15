%% script for bif and cont of spirals. OUTDATED, i.e., RW spirals should now rather be 
% computed (AND CONTINUED IN DELTA (dom-size) ) as relative equilibria. 
% Rather see cmds3.m, which also deals with mRWs, i.e., meandering spirals
close all; format compact; keep pphome;  
%% init, and find bifpoints from trivial branch, same HBPs as in cmds1
% but with al=1, thus recompute bifpoints to have the correct NF coeff
ndim=2; dir='tr-b'; p=[]; lx=pi; nx=30; par=[-0.25 0 1]; 
p=rotinit(p,nx,par); p=setfn(p,dir); p.nc.neig=50; 
p.sol.ds=0.01; p.nc.dsmax=0.01; p.nc.ilam=1; p.nc.bisecmax=20; p=cont(p,60); 
%% This may take some time; you might want to change bnr=1:7 to, e.g., bnr=2;
para=4; ds=0.4; dsmax=1.7; xi=1e-2; nsteps=50; lammax=3.1; 
flcheck=0; % we generally use flcheck=0. flcheck=1 only for s2 (most interesting) 
figure(2); clf; aux=[]; aux.tl=20; usrlam=[0 1 2 3];
for bnr=2; %[2 3 5 6 7] 
switch bnr
 case 1; p=hoswibra('tr-b','hpt1',ds,para,'s1',aux); 
 case 2; p=hoswibra('tr-b','hpt2',ds,para,'s2',aux); flcheck=1;
 case 3; p=hoswibra('tr-b','hpt3',ds,para,'s3',aux); flcheck=0; 
 case 4; p=hoswibra('tr','hpt4',ds,para,'s4',aux); 
 case 5; p=hoswibra('tr-b','hpt5',ds,para,'s5',aux); 
 case 6; p=hoswibra('tr-b','hpt6',ds,para,'s6',aux);  
 case 7; p=hoswibra('tr-b','hpt7',ds,para,'s7',aux);    
end 
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.nc.lammax=lammax; 
p.u(p.nu+2:p.nu+3)=[1;1]; % spiral case 
p.hopf.flcheck=flcheck; %p.file.smod=1; % switch on/off floquet (if somewhat slow)
AMG=1; p.sw.verb=3;  % set AMG=1 if ilupack available
if ~AMG; p=setbel(p,1,1e-3,10,@lss); % use BEL without ilupack 
else p=setilup(p,1e-3,200); %p.fuha.lss=@lss; p.fuha.blss=@lssAMG; 
    p=setbel(p,1,1e-3,10,@lssAMG);     
    p.usrlam=[]; % avoid problem with prec size
end 
%p.fuha.blss=@mbel; p.nc.mbw=2; p.hopf.ilss=0; % set to 0 if no ilupack 
p.usrlam=usrlam; tic; p=cont(p,nsteps);  toc
end 
%% continue in diff constant (domain size)
p=hoswiparf('s2','pt29','s2d',3,0.05); clf(2); 
p.usrlam=[0.1 0.2 0.3 0.4 0.5]; p.nc.lammin=0.09; p=cont(p,40);
%% BD of cont in delta
bpcmp=7; figure(3); clf; plotbra('s2d',3,bpcmp,'lab',[29,44],'cl','r'); 
xlabel('\delta'); ylabel('||u||_*'); 
%% continue in r again
p=hoswiparf('s2d','pt29','s2dr',1,-0.1); clf(2); p.sw.verb=2; 
p.usrlam=[0 1 2 3]; p.nc.lammin=-1; p=cont(p,40);
%% BD of cont in r
bpcmp=7; figure(3); clf; plotbraf('s2dr',3,bpcmp,'lab',[39 43],'cl','b'); 
xlabel('r'); ylabel('||u||_*'); 
%% BD of six branches, L^2
bpcmp=7; figure(3); clf; 
plotbra('tr',3,bpcmp,'cl','k'); plotbra('s1',3,bpcmp,'cl','b'); 
plotbra('s2',3,bpcmp,'lab',[12,15,29],'cl','r'); 
plotbra('s3',3,bpcmp,'cl','m'); plotbra('s4',3,bpcmp,'cl','b'); 
plotbra('s5',3,bpcmp,'cl','r'); plotbra('s6',3,bpcmp,'cl','m'); 
plotbra('s7',3,bpcmp,'cl','b'); 
axis([-0.3 3 0 1.2]); xlabel('r'); ylabel('||u||_*'); 
%% soln plot, here only plot profiles 
pstyle=12; proplotf('s2','pt12',1,1,pstyle); pause; proplotf('s2','pt29',1,1,pstyle); pause; 
proplotf('s3','pt34',1,1,pstyle); pause; proplotf('s4','pt12',1,1,pstyle); pause; 
proplotf('s5','pt20',1,1,pstyle); pause; proplotf('s6','pt18',1,1,pstyle); pause; 
proplotf('s7','pt18',1,1,pstyle);
%% profiles from cont in delta
proplotf('s2d','pt29',1,1,pstyle); pause; proplotf('s2d','pt44',1,1,pstyle); 
%% profiles for subsequent cont in r
proplotf('h2d-rc','pt39',1,1,pstyle); pause; proplotf('h2d-rc','pt43',1,1,pstyle); 
%% a-posteriori plots (for 3D plots) 
for i=1:4 
  switch i
    case 1; d='s2d'; pt='pt29'; 
    case 2; d='s2d'; pt='pt44'; 
    case 3; d='s2dr'; pt='pt39';    
    case 4; d='s2dr'; pt='pt43';   
  end
  p=loadp(d,pt); t1=['u_1 at ' d '/' pt]; t2=['u_2 at ' d '/' pt]; 
  u=p.hopf.y(1:p.nu,1); plotsolu(p,u,1,1,3); colormap cool; zlabel(''); 
  view(30,60); set(gca,'XTick',[]); set(gca,'yTick',[]); title(t1); 
  plotsolu(p,u,2,2,3); colormap cool; zlabel(''); view(30,60); 
  set(gca,'XTick',[]); set(gca,'yTick',[]); title(t2); pause
end
%% time integration, preparation: choose dir/pt, i.e., uncomment the desired one 
cmp=1; pstyle=2; wnr=2; notic=1; si=0; incr=100; nt=1000; 
lay=[2 2]; nplot=3; tsax=[0 50 0 4]; 
d='s3'; pt='pt34'; outdir='s3-34tint'; 
d='s5'; pt='pt18'; outdir='s5-18tint'; 
d='s2d'; pt='pt29'; outdir='s2d-29tint'; 
d='s2d'; pt='pt44'; outdir='s2d-44tint'; 
p=loadp(d,pt); hoplotrot(p,1,1,2); 
p.u(1:p.nu)=p.hopf.y(1:p.nu,1); u0=p.u(1:p.nu); p=setfn(p,outdir);
ts=[]; t0=0; npp=40; pmod=50; smod=5; tsmod=1; nc=0; 
p.mat.K=p.u(p.nu+3)*p.mat.K; % multiply K by delta for tint
%% call tint and plot time-series 
[p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); 
figure(4); clf; plot(ts(1,:), ts(2,:)); % point-vals
figure(5); clf; plot(ts(1,:), ts(3,:)); % difference in norm
set(gca,'FontSize',p.plot.fs); axis tight; 
xlabel('t'); ylabel('||u(t)-u_0||_{\infty}'); 
%% soln snapshots of time integration
ind=si:incr:si+nplot*incr;
tintplot2d(outdir,ind,wnr,cmp,pstyle,lay,notic); colormap cool;
figure(5); axis(tsax); 
%% Floquet a posteriori, choose points of interest! 
aux.nfloq=20; [muv1,~,~,~,h]=floqap('s2','pt15',aux); 