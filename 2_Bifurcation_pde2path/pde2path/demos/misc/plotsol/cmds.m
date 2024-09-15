close all; keep pphome;
%%
% creating 1D files
p1=[];p1=schnakinit(p1,1);p1.nc.nsteps=10000;p1=findbif(p1,1);
p1=swibra('h1','bpt1','s1');p1=cont(p1,20);
%% creating 2D files
p2=[];p2=schnakinit(p2,2);p2.nc.nsteps=10000;p2=findbif(p2,2);
p2=swibra('h2','bpt2','s2');p2=cont(p2,3);
%% creating 3D files
p3=[];p3=schnakinit(p3,3);p3.nc.nsteps=10000;p3=findbif(p3,1);
p3=swibra('h3','bpt1','s3');p3=cont(p3,2);