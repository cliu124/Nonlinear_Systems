close all; keep pphome; 
%% Schnakenberg cone, a=0.5 
dir='h2'; nx=50; ell=1.2; eps=0.05; lam=3.225; sig=0; D=60; a=0.5; s=0; dsmin=0.005; 
p=[]; par=[lam sig D a s eps]; p=schnakcinit(p,nx,par,ell); p.sol.ds=-dsmin; 
p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; 
% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%% hom.branch 
p=p0; p.nc.neig=20; p.sw.jac=1; p=cont(p,10);
%% BP1, use qswibra to compute kernel, then gentau 
bp='bpt1'; aux=[]; aux.m=4; aux.isotol=0.25; aux.besw=0; 
p0=qswibra('h2',bp,aux); p0.sw.bifcheck=0; p0.nc.dsmax=0.05; 
p0.pm.resfac=1e-2; p0.sol.ds=0.1; 
%% generating 3 branches via gentau 
p=gentau(p0,[1],'b1'); p=cont(p,5); coneplot(p); p=pmcont(p,20); 
p=gentau(p0,[0 1],'b2'); p=cont(p,5); coneplot(p); pause; p=pmcont(p,20); 
p=gentau(p0,[0 0 0 1],'b3'); p=cont(p,5); coneplot(p); pause; p=pmcont(p,20); 
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=7; plotbra('h1','pt20',fnr,pcmp,'cl','k','lsw',0);
plotbra('b1',fnr,pcmp,'cl','r','lab',10,'lp',15); 
plotbra('b2',fnr,pcmp,'cl','m','lab',10); 
plotbra('b3',fnr,pcmp,'cl',p2pc('b3'),'lab',10); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
axis([3 3.225 3 4.4]); 
%% solution plots 
p=loadp('b1','pt10'); plotsol(p); nola; coneplot(p,1); pause 
p=loadp('b2','pt10'); plotsol(p); nola; coneplot(p,1); pause 
p=loadp('b3','pt10'); plotsol(p); nola; coneplot(p,1); 
%% cont in a, b1
p=swiparf('b1','pt10','b1a',4); p.usrlam=[0.25 0.5]; p.sol.ds=0.05; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-6; p.sw.foldcheck=0; 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.1; p=pmcont(p,40); 
%% plot BD in a 
fnr=3; figure(fnr); clf; pcmp=0; 
plotbra('b1a','pt40',fnr,pcmp,'cl','r','lab',[0 30 40]); ylabel('||u_1||_2'); 
%% solution plots 
dir='b1a'; iset=[0 30 40]; %dir='b2a'; iset=[0 20]; %dir='b3a'; iset=[0 40];
for i=1:length(iset)
    pt=['pt' mat2str(iset(i))]; 
    p=loadp(dir,pt); plotsol(p); coneplot(p,2); pause 
end 
