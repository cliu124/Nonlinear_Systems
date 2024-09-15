%% 2D CPs at R=20; this is expensive (typical CP about 10Min) 
close all; keep pphome;
%% FSS->hex 
p=[]; p=ocinit(p,'2Df','pt29','hex','pt17');
p.oc.T=100; p.oc.nti=20; p.tomopt.Nmax=200; p.oc.verb=2; p.oc.msw=0;
%p.oc.sigmin=1e-5; p.oc.sig=1e-1; p.oc.sigmax=1; 
p.oc.tadevs=0.05; % allowed deviation from target CSS; 
alvin=[0.1 0.5 1]; tic; p1=isc(p,alvin); toc 
%% plot CP 
figure(15); clf; subp=2; ps=3; pfak=20; v=[20,60]; 
cpplot2D(p1,v,subp,ps,pfak,poc2fn(p)); 
%% FSS-> vertical stripes, like 1D 
p=[]; p=ocinit(p,'2Df','pt29','vs','pt22');
p.oc.T=100; p.oc.nti=20; p.tomopt.Nmax=200; p.oc.verb=2; p.oc.msw=0;
%p.oc.sigmin=1e-5; p.oc.sig=1e-1; p.oc.sigmax=1; 
p.oc.tadevs=0.05; alvin=[0.1 0.5 1]; tic; p2=isc(p,alvin); toc 
%% plot CP 
figure(15); clf; subp=2; ps=3; pfak=20; v=[20,60]; 
cpplot2D(p2,v,subp,ps,pfak,poc2fn(p)); 
%% FSS-> horizontal stripes, 
p=[]; p=ocinit(p,'2Df','pt29','hs','pt20');
p.oc.T=150; p.oc.nti=20; p.tomopt.Nmax=200; p.oc.verb=2; p.oc.msw=0;
%p.oc.sigmin=1e-5; p.oc.sig=1e-1; p.oc.sigmax=1; 
p.oc.tadevs=0.2; alvin=[0.1 0.5 1]; tic; p3=isc(p,alvin); toc; 
%% plot CP 
figure(15); clf; subp=2; ps=3; pfak=10; v=[20,60]; 
cpplot2D(p3,v,subp,ps,pfak,poc2fn(p)); 