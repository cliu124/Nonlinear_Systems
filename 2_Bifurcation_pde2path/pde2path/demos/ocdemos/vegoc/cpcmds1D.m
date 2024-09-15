close all; keep pphome;
%% R=20:
p=[]; p=ocinit(p,'FSS','pt6','p1','pt49'); p.oc.T=50; p.oc.nti=20; p.oc.msw=1; 
p.oc.mtom=0; % use hdv-bvp solver 
alvin=[0.1 0.25 0.5 0.75 1];v=[15,30]; p=isc(p,alvin);
figure(15); clf; vegsolplot(p,v,4,poc2fn(p));
%% R=10 
p=[]; p=ocinit(p,'FSS','pt22','p1','pt65'); p.oc.T=50; p.oc.nti=50; 
p.oc.msw=1; p.oc.mtom=0; alvin=[0.1 0.25 0.5 0.75 1]; p=isc(p,alvin);
figure(15); clf; v=[15,30]; vegsolplot(p,v,30,poc2fn(p));