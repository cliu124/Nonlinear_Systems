close all; keep pphome; 
%% Schnakenberg on pointed cone (a=0.2), scaled by 1/eps 
dir='h6'; nx=19; ell=1.2; eps=0.04; lam=3.225; sig=0; D=60; a=0.2; s=0; dsmin=0.005; 
p=[]; par=[lam sig D a s eps]; p=schnakcinit(p,nx,par,ell); p.sol.ds=-dsmin;  
p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.lammin=3; p.nc.usrlam=[3 3.1 3.2]; 
% dummy refine near 0 
p.nc.ngen=1; p.nc.maxt=20000; p.nref=100; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); 
p.tau=p.u; p0=p; 
%%
p=p0; p.t2sinterpol=1; p=tri2six(p); p.np, p.Kfn='K1.mat'; try delete(p.Kfn); catch; end; 
p.sw6=1; p.isw=1; p=oosetfemops(p); p.sol.xi=1/p.nu; plotsol(p,1,1,1); 
p0=p; 
%% hom.branch 
p=p0; p.nc.neig=20; p=findbif(p,2); p.sw.bifcheck=0; p=cont(p,10); 
%% first 2 branches 
p=swibra('h6','bpt1','c1',0.05); pause; p=cont(p,5); coneplot(p); p=pmcont(p,10); pause 
p=swibra('h6','bpt2','c2',0.02); pause; p=cont(p,5); coneplot(p); p=cont(p,20); 
%% 2ndary bifs 
p=swibra('c1','bpt1','c1-1',0.05); pause; p.nc.dsmax=0.05; p=cont(p,5); coneplot(p); p=pmcont(p,20); 
%%
p=swibra('c2','bpt3','c2-3',0.1); p=cont(p,5); coneplot(p); p=cont(p,20); 
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=7; plotbra('h6',fnr,pcmp,'cl','k','lsw',0);
plotbra('c1','pt20',fnr,pcmp,'cl','b','lab',15); 
plotbra('c2','pt30',fnr,pcmp,'cl','m','lab',20); 
plotbra('c1-1','pt20',fnr,pcmp,'cl',p2pc('o1'),'lab',20); 
plotbra('c2-3','pt20',fnr,pcmp,'cl',p2pc('r3'),'lab',20); 
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
