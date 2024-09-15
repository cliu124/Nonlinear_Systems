%% extension of bru1dcmds wrt bifs from periodic orbit 1dh1
close all; keep pphome; 
%% reload point from 1dh1 and run with more bisectionss
p=loadp('1dh1','pt10','1dh1b'); p.hopf.bisec=8; huclean(p); p=cont(p,20); 
%% bifurcations FROM 1dh1
ds=0.5; aux.sw=1; p=poswibra('1dh1b','bpt1','t1',ds,aux); 
p.sw.bifcheck=0; p.hopf.fltol=1e-2; % increase fl-tol due to large amplitude 
p.nc.tol=1e-3; p=cont(p,1); % do 1 step with large tol to get on bif.branch 
p.nc.tol=1e-8; p=cont(p,19); % decrease tol and continue further 
%% 2nd BP 
ds=0.5; aux.sw=1; p=poswibra('1dh1b','bpt2','t2',ds,aux); 
p.sw.bifcheck=0; p.hopf.fltol=1e-2; p.nc.tol=1e-4; p=cont(p,2); 
p.nc.tol=1e-8; p=cont(p,100); 
%% plot BD, L^2
bpcmp=9; wnr=4; figure(wnr); clf; plotbra('hom1d','pt45',wnr,bpcmp,'cl','k'); 
plotbra('1dh1','pt30',wnr,bpcmp,'cl',[1 0.5 0]); 
plotbra('t1','pt30',wnr,bpcmp,'cl','r','lab',[10],'fp',0); 
%plotbra('t2','pt50',wnr,bpcmp,'cl','m','lab',[70 80 90],'fp',0); 
plotbra('t2','pt50',wnr,bpcmp,'cl','m','lab',[20 30 40],'fp',0); 
xlabel('b'), ylabel('||u||_\infty'); %axis([2.75 3.25 2.75 5]); 
%% soln and floquet plotting 
aux.nfloq=100; v=[15,60];
dir='1dh1'; pt='bpt1'; hoplotf(dir,pt,1,1); figure(1); title(['u1 at ' dir '/' pt]); 
colormap cool; shading interp; view(v); muv1=floqap(dir,pt,aux); pause 
pt='bpt2'; hoplotf(dir,pt,1,1); figure(1); title(['u1 at ' dir '/' pt]); 
colormap cool; shading interp; view(v); muv1=floqap(dir,pt,aux); pause
dir='t1'; pt='pt10'; hoplotf(dir,pt,1,1); figure(1); title(['u1 at ' dir '/' pt]); 
colormap cool; shading interp; view(v); muv1=floqap(dir,pt,aux); pause
%%
dir='t2'; pt='pt20'; hoplotf(dir,pt,1,1); figure(1); title(['u1 at ' dir '/' pt]); 
colormap cool; shading interp; view(v); muv1=floqap(dir,pt,aux); pause
pt='pt30'; hoplotf(dir,pt,1,1); figure(1); title(['u1 at ' dir '/' pt]); 
colormap cool; shading interp; view(v); muv1=floqap(dir,pt,aux); pause
pt='pt40'; hoplotf(dir,pt,1,1); figure(1); title(['u1 at ' dir '/' pt]); 
colormap cool; shading interp; view(v); muv1=floqap(dir,pt,aux); 