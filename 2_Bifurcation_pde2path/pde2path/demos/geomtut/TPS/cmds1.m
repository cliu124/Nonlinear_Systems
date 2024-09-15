%% TPMS cont in scaling  
close all; keep pphome; 
%%
global p2pglob; p2pglob.vi=[20 15]; p2pglob.edc='k'; p2pglob.faceal=1; p2pglob.cut=0; 
p2pglob.tsw=0; p2pglob.showbd=0; p2pglob.pbctol=1e-6; p2pglob.cm='parula'; p2pglob.cb=0; 
lx=pi; ly=pi; lz=pi; ac=0.3; % triangulation width; other pars will be SET in Pinit
h0=0; v0=0; A0=0; sx=0; sy=0; sz=0; del=0; par=[h0; v0; A0; del; sx; sy; sz]; 
p=[]; p=Pinit(p,lx,ly,lz,ac,par); p=setfn(p,'P'); huclean(p); pplot(p); 
%% P to smaller delta 
p.file.smod=1; p=cont(p,15); 
%% Other direction
p=loadp('P','pt2','Pb'); p.nc.neig=20; p.sol.ds=-p.sol.ds; p=cont(p,20);
%% swibra 
p=swibra('P','bpt1','P1',-0.01); p.nc.tol=1e-6; pause; p=cont(p,20);
%% branch plot
f=3; c=[4 9]; mclf(f); plotbra('P','pt15',f,c,'cl','k','lab',[1 15]); 
plotbra('Pb','pt20',f,c,'cl',p2pc('gr1'),'lab',[7 20]);
plotbra('P1','pt20',f,c,'cl','m','lab',[7 20]); ylabel('A');
%% soln plot
plotsol('P','pt1'); pause; plotsol('P','pt15');pause; 
plotsol('Pb','pt7'); pause; plotsol('Pb','pt20'); pause 
plotsol('P1','pt7'); pause; plotsol('P1','pt20');
%% plotting some most unstable eigenvectors 
p2pglob.cb=1; p2pglob.vi=[20, 30]; p2pglob.showbd=2; 
ploteigsys('P','bpt1',6); pause; 
ploteigsys('Pb','pt7',6); pause; 
ploteigsys('P1','pt20',10); 
%% cont in H
p=loadp('P','pt1','PH'); p=resetc(p); p.sol.ds=0.1; mclf(2); 
p.nc.ilam=[1 5 6 7]; p.nc.nq=3; p.fuha.qf=@qf; p.fuha.qfder=@qjac; 
%%
p=cont(p,10); 
%% ------------ some tests:  refine, seems quite OK 
p2pglob.edc='k'; p2pglob.pbctol=1e-2; p.sw.orgper=0;
p=loadp('P','pt10','Pr'); 
plotsol(p); sig=0.1; p.sw.rlong=1; 
p=refineX(p,sig); p=cont(p,10); 
%% degcoarsen OK
p2pglob.pbctol=1e-4; p=loadp('Pb','pt15','Pbc'); pplot(p,1); nit=2; keepbd=1; 
sig=0.1; p=degcoarsenX(p,sig,nit,keepbd); p.sol.ds=-p.sol.ds; pause; p=cont(p,5); 
%% repeat degcoarsen, OK 
p2pglob.pbctol=1e-10; p=loadp('Pbc','pt20','Pbc'); sig=0.05; 
p=degcoarsenX(p,sig,nit,keepbd); pause; p=cont(p,10); 