%% SLOC canonical paths. CP from p3/pt19 to f1/pt13; nat.continuation
p=[]; p=ocinit(p,'p3','pt19','f1','pt13'); % set standard options, also set
% p3/pt19 as start and f1/pt13 as end point of the canonical path
p.oc.rhoi=1; % index of discount rate rho, rho is at p.u(p.nu+p.nq+rhoi)
p.oc.T=100; % truncation time
% reset some oc-params to customized values 
p.oc.nti=41; % initial # of points in t-mesh 
p.oc.msw=1; % use secant predictors (after first step) in isc
alvin=[0.25 0.5 0.75 1]; % desired alpha-values; these can also be split, 
% e.g., first call isc with alvin=[0.25 0.5], then again with alvin=[0.75 1]  
p=isc(p,alvin); va=[15,30]; slsolplot(p,va); % continuation call and plot 
%% start with small T, then set T free (test T adaptation) 
p=[]; p=ocinit(p,'p3','pt19','f1','pt13'); p.oc.rhoi=1; p.oc.T=30; 
p.oc.nti=41; p.oc.msw=1; p.oc.tadevs=1e-2; p.oc.verb=2; p.oc.mtom=0; 
alvin=[0.25 0.5 0.75 1]; p=isc(p,alvin); 
va=[15,30]; slsolplot(p,va); tadev(p); % some plots 
%% decrease ||u(1)-\uhat|| (increasing T from 43 to about 50) 
p2=p; p2.oc.tadevs=1e-4; p2=isc(p2,1); slsolplot(p2,va); tadev(p2); 
%% 20-11-16: see cpdemo1D_leg for fold in alpha (old syntax)! 
%% arclength continuation, initial call
p=[]; p=ocinit(p,'p3','pt19','f1','pt13'); % gives al=1 in just 2 arclength steps 
%p=ocinit(p,'p1','pt15','f1','pt13'); % needs more arclength steps 
%p=ocinit(p,'p1','pt68','f1','pt13'); % for al-fold, see cpdemo1D_legacy variant 
p.oc.rhoi=1; p.oc.nti=41; p.oc.T=50; p.oc.mtom=0; 
p.oc.verb=2; p.oc.msw=1; p.oc.almax=0.2; p.tomopt.tol=1; % no it in 1st step, cause that lets al jump! 
p.oc.sig=0.05;  p.oc.sigmax=1; % stepsize and max stepsize 
p.tomopt.Stats_step='on'; % give mtom diagnostics 
p.oc.retsw=1; % return full alpha and soln history 
alvin=[0.01 0.1]; arcsteps=1; p=isc(p,alvin,arcsteps); % continuation call, startup 
%% subsequent arclength-continuation-calls (repeat this cell to continue further) 
p.tomopt.tol=1e-8;
p=isc(p,[],20); 
%% plots, including a simple plot of J over alpha, 
figure(6); clf; plot(p.hist.alpha(1,2:end),p.hist.vv(1,2:end),'-*'); 
xlabel('\alpha'); ylabel('J_{a}'); set(gca,'FontSize',14); 
va=[25,15]; slsolplot(p,va); 
%% fix al from arc cont to some given value and compute CP
os=4; p.cp.u=p.hist.u{end-os}; p.cp.t=p.hist.tt{end-os};
p.cp.par=p.hist.par{end-os}; al=p.cp.par(2);% extract sol from history
al=round(al*10)/10
p=isc(p,al); v=[100,30]; slsolplot(p,va); % compute at desired al value and plot 
% here gives different J_c than the arclength solution! need to be fixed