%% c1: some mesh-refinement, here just for one fixed solution
p=loadp('b1','pt10','b1'); p.nc.maxt=100; 
%p.fuha.e2rs=@e2rs_ad_hoc; % uncomment to use ad hoc element selector
p=meshada(p,'ngen',3,'sig',0.5); 
%% c2: continuation with mesh-adaption depending on error bound 
p=swibra('tr','bpt1','b1ref',0.1); p.nc.sig=0.95; % large sig gives fine mesh! 
p.nc.lammax=2; p.sw.errcheck=2; p.nc.errbound=0.2; p=cont(p,40);  
figure(3); clf; plotbra(p,3,-1,'lsw',0,'labi',5); % plot error-est on branch
%% c3: continuation with mesh-adaption each amod-th step 
p=swibra('tr','bpt1','b1ref',0.1); p.nc.sig=0.95; 
p.sw.errcheck=0; p.nc.amod=5; p.nc.maxt=200; p.nc.lammax=2;  p=cont(p,40);  
figure(3); clf; plotbra(p,3,-1,'lsw',0,'labi',3); % plot error-est on branch
plotsol(p,1,1,1,'pstyle','*');
%% c3b: legacy setup of meshadac, flagged by p.sw.scoarse=0; 
p=swibra('tr','bpt1','b1ref',0.1); p.nc.sig=0.95; p.sw.scoarse=0; 
p.sw.errcheck=0; p.nc.amod=5; p.nc.maxt=200; p.nc.lammax=2;  p=cont(p,40);  
figure(4); clf; plotbra(p,4,-1,'lsw',0,'labi',3); % error-est=penultimate entry of bradat
plotsol(p,1,1,1,'pstyle','*');