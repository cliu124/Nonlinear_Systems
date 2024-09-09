%% commands for Schnakenberg on longer domain, with snaking (not yet on octave) 
close all; keep pphome; 
p=[]; kc=sqrt(sqrt(2)-1); lx=8*pi/kc; ly=lx/sqrt(3)/4; nx=100; par=[3.3, 0, 60]; 
sw.sym=2; p=schnakinit(p,[lx,ly],nx,par,sw); p.np % init with criss-cross mesh 
p.fuha.outfu=@schnakbra; p=setfn(p,'bhom'); p=cont(p,10);
%% hex+ via qswibra
aux=[]; aux.m=4; p0=qswibra('bhom','bpt1',aux); p0.sol.ds=-0.1; p0.sw.bifcheck=1;
p=seltau(p0,1,'bh+',2);  p=pmcont(p,30); 
%% stripes+ via gentau 
p=gentau(p0,[0 0 0 1]); p=setfn(p,'bs+'); p=cont(p,20);
%% plot BD 
fnr=3; figure(fnr); clf; cmp=6; 
plotbra('bhom','bpt2',fnr,cmp,'cl','k'); plotbra('bs+','pt20',fnr,cmp,'cl','b'); 
plotbra('bh+','pt20',fnr,cmp, 'cl','m'); 
%% soln plot 
p=loadp('bs+','pt20'); plotsol(p); colormap hot; fouplot(p,10,1,[20,3],'sn1/pt20'); 