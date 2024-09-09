%% SWs and RWs for GKS2000 RD model for spirals; slightly outdated, i.e., RW are now rather 
% computed as relative equilibria; 
close all; format compact; keep pphome; 
%% C1: init, and find bifpoints from trivial branch
ndim=2; dir='tr'; p=[]; lx=pi; nx=30; par=[-0.22 0 1 0]; 
p=rotinit(p,nx,par); p=setfn(p,dir); p.nc.neig=50; p.sol.ds=0.01; p.nc.dsmax=0.01; 
p.nc.ilam=1; p.nc.bisecmax=10; p=cont(p,30); 
%% C2: mode plot at bif 
hospatplot('tr','hpt2',6,1,3);
%% C3: first 7 branches; at HBPs with angular wnr 0 compute SWs, otherwise RWs
% this may take some time; you might want to change bnr=1:7 to, e.g., bnr=1;
% RWs are rigid rotations, hence can also be computed via twswibra, see next cell
% SWs can be more eff. be computed with a phase condition, see, e.g., cgldisk 
para=4; ds=0.4; dsmax=1; xi=1e-2; nsteps=10; 
figure(2); clf; aux=[]; aux.tl=20; aux.pstyle=3; 
hoaux.lay=[1 4]; hoaux.pind=[1 4 7 10]; 
hoaux.xtics=''; hoaux.ytics=''; hoaux.view=[0 90]; aux.hoaux=hoaux;
for bnr=[1 4]
 switch bnr
 case 1; p=hoswibra('tr','hpt1',ds,para,'sw1',aux); 
 case 2; p=hoswibra('tr','hpt2',ds,para,'rw2',aux); 
 case 3; aux.z=[1 1i];  p=hoswibra('tr','hpt3',ds,para,'rw3',aux);  % needs some tweaking 
 case 4; aux.z=1; p=hoswibra('tr','hpt4',ds,para,'sw4',aux); 
 case 5; aux.z=[1 1i]; p=hoswibra('tr','hpt5',ds,para,'rw5',aux); 
 case 6; aux.z=1; p=hoswibra('tr','hpt6',ds,para,'rw6',aux); 
 case 7; p=hoswibra('tr','hpt7',ds,para,'rw7',aux);  
 end 
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.usrlam=[]; p.nc.tol=1e-6; 
p.u(p.nu+2:p.nu+3)=[0;1]; % default case 
p.hopf.flcheck=0; p.file.smod=1; % switch off floquet (slow)
AMG=1; p.sw.verb=3;  % set AMG=1 if ilupack available
if ~AMG; p=setbel(p,1,1e-3,10,@lss); % use BEL without ilupack 
else p=setilup(p,1e-3,100); p.fuha.lss=@lss; p.fuha.blss=@lssAMG; 
     % no bel on this one, cause ilupack seems indifferent to borders! 
end 
tic; p=cont(p,nsteps);  toc
end
%% C4: compute RWs via twswibra
aux=[]; spar=4; % speed (here rot) parameter index 
for bnr=[5 6]; %[2 3 5 6 7]; 
switch bnr
 case 2; p=twswibra('tr','hpt2',spar,1,'rw2b',aux); 
 case 3; p=twswibra('tr','hpt3',spar,2,'rw3b',aux); 
 case 5; p=twswibra('tr','hpt5',spar,3,'rw5b',aux);
 case 6; p=twswibra('tr','hpt6',spar,1,'rw6b',aux);
 case 7; p=twswibra('tr','hpt7',spar,4,'rw7b',aux);
end 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; 
p.u0x=p.mat.Krot*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.05*p.tau(1:p.nu); 
plotsolu(p,p.u0,6,1,3); plotsolu(p,5*p.u0x,11,1,3);
p.nc.nq=1; p.nc.ilam=[1;spar];  % 1 phase-cond, speed as second parameter
p.sw.bifcheck=0; p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; p.sw.verb=2; 
p.sw.bprint=4; clf(2); p.nc.dsmax=0.1; p.sol.ds=0.04; pause; p=cont(p,10); 
end
%% soln plot, here use specialized 'hoplot' ! 
pstyle=2;
hoplotrotf('sw1','pt10',4,1,pstyle); pause; hoplotrotf('rw2','pt10',4,1,pstyle); pause; 
hoplotrotf('rw3','pt10',4,1,pstyle); pause; hoplotrotf('sw4','pt10',4,1,pstyle); pause; 
hoplotrotf('rw5','pt10',4,1,pstyle); pause; hoplotrotf('rw6','pt10',4,1,pstyle); 
%% tint, preparations 
p=loadp('rw3','pt10'); hoplotrot(p,1,1,2); dir='tinth3'; 
p.u(1:p.nu)=p.hopf.y(1:p.nu,1); u0=p.u(1:p.nu); p=setfn(p,dir);
ts=[]; t0=0; npp=40; pmod=50; smod=5; tsmod=1; nc=0; 
%% call tint and plot time-series 
nt=200; [p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); 
figure(4); clf; plot(ts(1,:), ts(2,:)); % point-vals
figure(5); clf; plot(ts(1,:), ts(3,:)); % difference in norm
set(gca,'FontSize',p.plot.fs); axis tight; 
xlabel('t'); ylabel('||u(t)-u_0||_{\infty}'); 
% soln snapshots of time integration
dir='tinth3'; si=0; incr=70; nt=1*npp/smod; wnr=2; cmp=2; pstyle=2; 
ind=si:incr:11*incr; lay=[2 6]; notic=1; 
tintplot2d(dir,ind,wnr,cmp,pstyle,lay,notic);
%% Floquet a posteriori, choose points of interest! 
aux.nfloq=20; [muv1,~,~,~,h]=floqap('rw2','pt10',aux); 