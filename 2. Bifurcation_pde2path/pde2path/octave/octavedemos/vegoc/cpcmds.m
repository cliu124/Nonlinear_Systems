close all; keep pphome;
%% R=20: CP from FSS/pt29 to p1/pt49, set some options, then call iscnat and plot
p=[]; p=ocinit(p,'FSS','pt29','p1','pt49'); p.oc.T=50; p.oc.nti=20; p.oc.msw=1; 
p.oc.mtom=0; % use HdW solver, mtom has problems under octave! 
alvin=[0.1 0.25 0.5 0.75 1];v=[15,30]; p=isc(p,alvin);
figure(15); clf; vegsolplot(p,v,4,poc2fn(p));
%% R=10 
p=[]; p=ocinit(p,'FSS','pt45','p1','pt65'); p.oc.T=100; p.oc.nti=20; p.oc.mtom=0; 
p.oc.msw=1; alvin=[0.1 0.25 0.5 0.75 1]; p=isc(p,alvin);
figure(15); clf; v=[15,30]; vegsolplot(p,v,30,poc2fn(p));