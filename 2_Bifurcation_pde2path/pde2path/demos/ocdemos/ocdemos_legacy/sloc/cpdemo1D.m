% script for SLOC canonical paths (short version, see cpdemo1D_full for more)
global s0 s1 u0 u1 Psi par xi um1 um2 sig;  % set global vars 
%% Cell 1: preparations; put filenames into fn, set some bvp parameters 
sd0='f1'; sp0='pt13'; sd1='p3'; sp1='pt19'; flip=1; % p3->FSC
fn=setfnflip(sd0,sp0,sd1,sp1,flip); % filename-struct 
opt=[]; opt=ocstanopt(opt); % set options and params  for CP computations
% set oc-params for which there is no default: 
opt.rhoi=1; % index of discount rate rho, rho is at p.u(p.nu+p.nq+rhoi)
opt.t1=100; % truncation time
% reset some oc-params to customized values 
opt.start=1; % 1 for first call
opt.tv=[]; sol=[]; % t-mesh, solution, empty at startup 
opt.nti=20; % initial # of points in t-mesh
opt.retsw=0; % only return final step in iscnat 
opt.msw=1; % use secant predictors (after first step) in iscnat 
usec=[]; % first secant (naturally) empty 

alvin=[0.25 0.5 0.7 1]; % desired alpha-values; these can also be split, e.g., 
% first call iscnat with alvin=[0.25 0.5], then again with alvin=[0.7 1] 
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,sol,usec,opt,fn); % continuation call, 
% see iscnat source and tutorial for more explanation 
va=[15,30]; % viewing angles used in solution plots
slsolplot(sol,va); % some plots 
%% Cell 2: a fold in alpha, here iscarc needed; Prep. and initial iscarc call
sd0='f1'; sp0='pt13'; sd1='p1'; sp1='pt68'; flip=1; usec=[]; 
fn=setfnflip(sd0,sp0,sd1,sp1,flip); 
esol=[]; % no extended sol at startup (iscarc on startup calls isnat anyway) 
opt.start=1; % 1 for first call 
opt.alvin=[0.2 0.3]; % alpha values for startup of iscarc 
opt.nsteps=3; % # of steps for iscarc 
sig=0.1;  opt.sigmax=1; % stepsize and max stepsize 
opt.Stats_step='on'; % give mtom diagnostics 
opt.retsw=1; % return full alpha and soln history 
[alv,vv,usec1,esol1,tlv,tv,uv]=iscarc(esol,usec,opt,fn); % continuation call 
opt.start=0;  % set to 0 for subsequent calls 
%% Cell 3: subsequent iccarc-calls (repeat this cell to continue further) 
opt.nsteps=35; usec=usec1; esol=esol1; % new input (for repeated calls) 
[alv1,vv1,usec1,esol1,tlv1,tv1,uv1]=iscarc(esol,usec,opt,fn); 
% append output to previous steps for repeated calls 
alv=[alv alv1]; vv=[vv vv1]; tlv=[tlv tlv1]; tv=[tv; tv1]; uv=[uv; uv1]; 
alv0=alv; vv0=vv; tv0=tv; uv0=uv; tlv0=tlv; %  save results for skibademo.m 
%% Cell 4: a simple plot of J over alpha  
figure(6); clf; plot(alv(1,:),vv(1,:),'-*'); set(gca,'FontSize',14); 
xlabel('\alpha','FontSize',14); ylabel('J_{a}','FontSize',14);
%% Cell 5: plot selected (here: the last) path from continuation
j=length(tlv); tl=tlv(j); va=[25,15]; 
n=s1.nu; sol.x=tv(j,1:tlv(j));sol.y=squeeze(uv(j,1:n,1:tlv(j)));
al=alv(j), u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); u1=s1.u(1:s1.nu); slsolplot(sol,va); 
%% Cell 6: fix al from iscarc to some given value and compute CPs
j=22; alv(j), n=s1.nu; sol.x=tv(j,1:tlv(j));sol.y=squeeze(uv(j,1:n,1:tlv(j))); 
al=0.6; u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); u1=s1.u(1:s1.nu); 
opt.M=s1.mat.M; sol=mtom(@mrhs,@cbcf,sol,opt); v=[100,30]; slsolplot(sol,va); 