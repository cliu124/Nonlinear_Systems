close all; format compact; keep pphome; % clean up 
%% init, and continuation of trivial branch 
p=[]; lx=pi; nx=30; p=cGLinit(p,lx,nx); % initialize  
% parameters: 1=r,2=speed,3=nu,4=mu,5=c3,6=c5,7=gam,8=sym,9=domain scale
par=[-1; 0; 0;  0;  -1; 1; 0;   0;   1];  
%     r  s  nu  mu  c3 c5  gam  sym  l
u=zeros(p.np,1); v=u; p.u=[u;v; par]; % initial guess (here trivial) and pars 
dir='zero'; p=setfn(p,dir); p.fuha.savefu(p); % set dirname and save
%% PERIODIC BC: continue zero solution on periodic domain
p=loadp('zero','pt0','per/zero'); p=box2per(p,1); % switch to periodic BC
p.sw.eigsstart=1; % set v0\equiv 1 for inv.vector iteration (default).  
% Set p.sw.eigsstart=0 for random v0; this is useful here to get the full
% spectrum. Backdraw: uncontrolled kernel vectors and thus also branches 
p=cont(p,20); % continuation of (here) trivial branch, incl. bif-detec
%% switch to branch for a few steps without phase-cond (PC) 
p=swibra('per/zero','bpt2','per/ini',0.025); 
p.sw.bifcheck=0; p.sw.spcalc=0; % switch off bif-detection (is random without PC) 
p.file.smod=1; p=cont(p,4); % do 1 initial step (for check), and store all 
%% add phase conditions to 1st point of cell 3 and continue
p=loadp('per/ini','pt1','per/stand'); p=resetc(p); p.file.smod=10; 
p.nc.ilam=[1;2;3];  p.nc.nq=2; % 2 phase-cond, speed and nu as add parameters
p.fuha.qf=@qf_sym;      % function handle for aux. eqn.
p.sw.qjac=1; p.fuha.qfder=@qfder_sym; % analytical jac for aux. eqn.
p.sw.bprint=2; clf(2); p.nc.dsmax=0.05; p=cont(p,30);
%% add only transl. PC to 1st point of cell 3 -> no conv. or uncontr. rotations 
p=loadp('per/ini','pt1','per/rot0'); p=resetc(p); p.file.smod=10; 
p.nc.nq=1; p.nc.ilam=[1;2];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf_trans; p.sw.qjac=1; p.fuha.qfder=@qfder_trans; 
p.sw.bprint=2; clf(2); p.nc.dsmax=0.05; p.sol.ds=-0.05; p=cont(p,40);
%% point 5 of cell 3 (TW), only transl phase-cond sufficient 
p=loadp('per/ini','pt5','per/wtstand'); p=resetc(p); p.file.smod=10; 
p.nc.nq=1; p.nc.ilam=[1;2];  % 1 phase-cond, speed as second parameter
p.branch=[bradat(p); p.fuha.outfu(p,p.u)]; figure(2); clf; 
p.fuha.qf=@qf_trans; p.sw.qjac=1; p.fuha.qfder=@qfder_trans; 
p.sw.bprint=2; p.nc.dsmax=0.1; p.sol.ds=-0.1; p=cont(p,20);
%% Plotting cmds in plotcmds.m 