close all; keep pphome; 
%% commands for Schnakenberg on a small 2D domain with hex lattice
p=[]; kc=sqrt(sqrt(2)-1); lx=2*pi/kc; ly=lx/sqrt(3); par=[3.3, 0, 60]; 
nx=80;sw.sym=2; p=schnakinit(p,[lx,ly],nx,par,sw); p.np % init with criss-cross mesh 
p.pm.resfac=1e-4; p.sol.ds=-0.1; p=setfn(p,'hom'); p=cont(p,30);
%% hex via qswibra, continue in both directions
p0=qswibra('hom','bpt1'); p0.nc.dsmin=0.1; p0.sw.bifcheck=1; pause
p=seltau(p0,3,'h-',2); p.sol.ds=0.1; p=pmcont(p,20);
p=seltau(p0,3,'h+',2); p.sol.ds=-0.1; p=pmcont(p,20);
%% s+ and s- via gentau
p=gentau(p0,[0 1],'s-'); p.sol.ds=-0.05; p=pmcont(p,20); 
p=gentau(p0,[0 1],'s+'); p=pmcont(p,50); 
%% b+ and b-
p=swibra('s+','bpt1','b+',0.05); p.nc.dsmin=0.05; p.pm.resfac=1e-3; p.file.smod=5; p=pmcont(p,20);
p=swibra('s-','bpt1','b-',0.05); p.nc.dsmin=0.05; p.pm.resfac=1e-3; p.file.smod=5; p=pmcont(p,10);
%% plot BD 
fnr=3; figure(fnr); clf; plotbraf('hom','pt20',fnr,4,'cl','k');
plotbra('s+','pt40',fnr,4,'cl','b','lab',35); plotbra('s-','pt20',fnr,5,'cl','b','lab',35); 
plotbra('b+','pt20',fnr,4,'cl','r','lab',5,'lp',18); 
plotbra('b-','pt10',fnr,5,'cl','g','lab',5,'lp',8); 
plotbra('h+','pt20',fnr,4, 'cl','m','lab',20);  plotbra('h-','pt20',fnr,5, 'cl','m','lab',20); 
text(2.2,2,'hom','Fontsize',16,'color','k'); axis([2 3.25 0.7 6.2]); xlabel('\lambda'); ylabel('max/min(u_1)'); 
%% plot sol. incl. Fourier
p=loadp('b-','pt5'); plotsol(p); fouplot(p,10,1,[5 4],'b-/pt5'); pause
p=loadp('b+','pt5'); plotsol(p); fouplot(p,10,1,[5 4],'b+/pt5'); pause
p=loadp('h-','pt20'); plotsol(p); fouplot(p,10,1,[5 4],'h-/pt20'); pause
p=loadp('h+','pt20'); plotsol(p); fouplot(p,10,1,[5 4],'h+/pt20'); pause
plotsol('s+','pt20',1,1,2); pause; plotsol('s-','pt20',1,1,2); 