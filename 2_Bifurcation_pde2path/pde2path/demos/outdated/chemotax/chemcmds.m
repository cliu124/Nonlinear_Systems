%% demo chemotax
% quasilinear system, here run with numerical G_u
%% Init
close all; % clear all; 
p=[]; p=cheminit(p); screenlayout(p);  
p.plot.pstyle=2; p.plot.pcmp=1; p=cont(p);
%% 5 bif. branches, low tolerance, switching off bifdetec for speed 
p=swibra('p','bpt1','q1',0.5); p.nc.tol=1e-6; 
p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p); 
p=swibra('p','bpt2','q2',0.5); p.nc.tol=1e-6; 
p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p); 
p=swibra('p','bpt3','q3',0.5); p.nc.tol=1e-6; 
p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p); 
p=swibra('p','bpt4','q4',0.5); p.nc.tol=1e-6; 
p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p); 
p=swibra('p','bpt5','q5',0.5); p.nc.tol=1e-6; 
p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p); 
%% plot bifurcation diagram
figure(4);clf(4);cmp=1;
plotbraf('p',4,cmp,'cl','k');
plotbraf('q1',4,cmp,'cl','b');
plotbraf('q2',4,cmp,'cl','r');
plotbraf('q3',4,cmp,'cl','m');
plotbraf('q4',4,cmp,'cl','c');
plotbraf('q5',4,cmp,'cl','g');
axis([10 22 0 0.8]);xlabel('\lambda');
ylabel('||u_1-1||_{L^1}/|\Omega|');
%% plot some solns 
plotsolf('q1','pt10',1,1,2); 
plotsolf('q2','pt10',2,1,2); 
