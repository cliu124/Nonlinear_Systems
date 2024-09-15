close all; keep pphome; 
%% Schnakenberg on pointed cone (a=0.2), scaled by 1/eps 
dir='h1'; nx=50; ell=1.2; eps=0.04; lam=3.225; sig=0; D=60; a=0.2; s=0; dsmin=0.005; nx=20; 
p=[]; par=[lam sig D a s eps]; p=schnakcinit(p,nx,par,ell); p.sol.ds=-dsmin;  
p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.lammin=3; p.nc.usrlam=[3 3.1 3.2]; p0=p; 
% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%% hom.branch 
p=p0; p.nc.neig=20; p=findbif(p,3); %p.sw.bifcheck=0; p=cont(p,20); 
%% first 2 branches 
p=swibra('h1','bpt1','c1',0.05); p=cont(p,5); coneplot(p); p=pmcont(p,30); 
p=swibra('h1','bpt2','c2',0.02); p=cont(p,5); coneplot(p); p=pmcont(p,30); 
%% 2ndary bifs 
p=swibra('c1','bpt2','c1-1',0.1); p=cont(p,5); coneplot(p); p=pmcont(p,20); 
p=swibra('c1','bpt3','c1-2',0.1); p=cont(p,5); coneplot(p); p=pmcont(p,20); 
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=7; plotbra('h1','pt20',fnr,pcmp,'cl','k','lsw',0);
plotbra('c1',fnr,pcmp,'cl','b','lab',15); 
plotbra('c2','pt20',fnr,pcmp,'cl','m','lab',20); 
plotbra('c1-1','pt20',fnr,pcmp,'cl',p2pc('o1'),'lab',20); 
plotbra('c1-2','pt20',fnr,pcmp,'cl',p2pc('r3'),'lab',20); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
%% solution plots 
v=[20 70]; p=loadp('c1','pt15'); plotsol(p); nola; coneplot(p,1); view(v); pause 
p=loadp('c1-1','pt20'); plotsol(p); nola; coneplot(p,1);  view(v);pause 
p=loadp('c1-2','pt20'); plotsol(p); nola; coneplot(p,1);  view(v);pause 
p=loadp('c2','pt20'); plotsol(p); nola; coneplot(p,1);  view(v);
%% cont in a, often closed loops 
p=swiparf('c1-1','pt10','c1a',4); p.usrlam=[0.25 0.5]; p.sol.ds=0.05; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-8; p.sw.foldcheck=0; p.nc.intol=-1e-4; 
p.sw.bifcheck=0; p.nc.lammax=2; p.nc.lammin=0.1; p.nc.dsmax=0.1; p=cont(p,40); 
%% BD in a
fnr=3; figure(fnr); clf; pcmp=0; 
plotbra('c1a','pt50',fnr,pcmp,'cl','r','lab',[0 30 40]); ylabel('||u_1||_2'); 
%% solution plots 
dir='c1a'; iset=[0 30 40 50]; %dir='c2a'; iset=[0 20]; 
for i=1:length(iset)
    pt=['pt' mat2str(iset(i))]; 
    p=loadp(dir,pt); plotsol(p); coneplot(p,2); pause 
end 
