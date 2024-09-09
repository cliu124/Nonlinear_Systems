%% ac2D, test sfem=1 vs legacy impl. 
close all; keep pphome; 
%% init and find 1st BP with diferent nx 
lx=2*pi; ly=pi; p=[]; par=[1 -0.2 1 0]; % parameters [c lambda gamma d]
p=acinit(p,lx,ly,20,par); p=setfn(p,'tr20'); p.file.smod=0; p=findbif(p,1);
p=acinit(p,lx,ly,50,par); p=setfn(p,'tr50'); p.file.smod=0; p=findbif(p,1);
p=acinit(p,lx,ly,80,par); p=setfn(p,'tr80'); p.file.smod=0; p=findbif(p,1);
%% 
p=swibra('tr20','bpt1','t1a',0.1); p=cont(p); 
p=swibra('tr20','bpt1','t1b',0.05); p.cl=1; p=cont(p); 
%% 
p=swibra('tr50','bpt1','t2a',0.1); p=cont(p); 
p=swibra('tr50','bpt1','t2b',0.1); p.cl=1; p=cont(p); 
%% 
p=swibra('tr80','bpt1','t3a',0.1); p=cont(p); 
p=swibra('tr80','bpt1','t3b',0.1); p.cl=1; p=cont(p); 
%% check differences 
p1=loadp('t1a','pt17'); p2=loadp('t1b','pt17'); 
plotsolu(p1,p1.u-p2.u,1,1,1); getaux(p1)', getaux(p2)'
title('n_x=20'); 
%% check differences 
p1=loadp('t2a','pt18'); p2=loadp('t2b','pt18'); 
plotsolu(p1,p2.u-p1.u,1,1,1); getaux(p1)', getaux(p2)'
title('n_x=50'); 