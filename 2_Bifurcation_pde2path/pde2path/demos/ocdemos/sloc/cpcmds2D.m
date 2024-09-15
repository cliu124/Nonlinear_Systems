%% SLOC canonical paths. CP from h2/pt15 to 2DFSC/pt13; nat.continuation
p=[]; p=ocinit(p,'h2','pt15','2DFSC','pt13'); % set standard options, also set
% p3/pt19 as start and f1/pt13 as end point of the canonical path
p.oc.rhoi=1; % index of discount rate rho, rho is at p.u(p.nu+p.nq+rhoi)
p.oc.T=100; % truncation time
% reset some oc-params to customized values 
p.oc.nti=41; % initial # of points in t-mesh 
p.oc.msw=1; % use secant predictors (after first step) in isc
alvin=[0.25 0.5 0.75 1]; % desired alpha-values; these can also be split, 
% e.g., first call isc with alvin=[0.25 0.5], then again with alvin=[0.75 1]  
p=isc(p,alvin); % go 
%%
va=[15,30]; nsub=2; pstyle=3; cpplot2D(p,va,nsub,pstyle);  % plot CP  
%% start with small T, then set T free (test T adaptation) 
p=[]; p=ocinit(p,'h2','pt15','2DFSC','pt13'); p.oc.rhoi=1; p.oc.T=30; 
p.oc.nti=41; p.oc.msw=1; p.oc.tadevs=1e-2; p.oc.verb=2; 
alvin=[0.25 0.5 0.75 1]; p=isc(p,alvin);  tadev(p);
cpplot2D(p,va,nsub,pstyle);  % plot CP  
%% decrease ||u(1)-\uhat|| (increasing T from 39 to 48) 
p2=p; p2.oc.tadevs=1e-4; p2=isc(p2,1); 
cpplot2D(p2,va,nsub,pstyle); tadev(p2); 
