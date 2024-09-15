%% Schwarz P, cont in H 
close all; keep pphome; 
%%
global p2pglob; p2pglob.vi=[20 15]; p2pglob.edc='k'; p2pglob.faceal=1; p2pglob.cut=0; 
p2pglob.tsw=2; p2pglob.showbd=0; p2pglob.pbctol=1e-6; p2pglob.cm='parula'; p2pglob.cb=0; 
lx=pi; ly=pi; lz=pi; ac=0.2; % triangulation width; other pars will be SET in Pinit
h0=0; v0=0; A0=0; sx=0; sy=0; sz=0; del=0; par=[h0; v0; A0; del; sx; sy; sz]; 
p=[]; p=Pinit(p,lx,ly,lz,ac,par); p=setfn(p,'PH'); huclean(p); pplot(p); p.np
%%
p.nc.ilam=[1 5 6 7]; p.nc.nq=3; p.fuha.qf=@qf; p.fuha.qfder=@qjac; 
p.nc.mu2=0.04; mclf(2); p.sol.ds=-0.01; p.nc.eigref=-1; p.nc.neig=10; p.u(p.nu+1)=-1e-6; 
p.usrlam=0; p=cont(p,15); 
%% Other direction
p=loadp('PH','pt0','PHb'); p.sol.ds=-p.sol.ds; p=cont(p,15);
%% swibra
p=swibra('PH','bpt1','bPH1',-0.01); p.nc.tol=1e-8; pause; p=cont(p,20);
%% qswibra; ABE work, but looks like just 2 orientations come back 
aux.besw=1; aux.m=2; aux.ral=0; aux.isotol=1e-2; p1=qswibra('PH','bpt1',aux);
%% z-axis orient 
p=gentau(p1,[0.2 1],'za'); pause; p.sol.ds=-0.05; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,6);
p=gentau(p1,[0.2 1],'zb'); pause; p.sol.ds=0.05; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,15);
%% a possible extra branch without D4 symmetry? -using that predictor we fall back to D4! 
p=seltau(p1,2,'zm',2); pause; p.sol.ds=0.05; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,6);
%% swibra at positive H:
aux.besw=0; aux.m=2; aux.ral=1; aux.isotol=1e-4; p1=qswibra('PHb','bpt1',aux);
%% z-axis 
p=gentau(p1,[-0.5 1],'za2'); pause; p.sol.ds=-0.05; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,8);
p=gentau(p1,[-0.5 1],'zb2'); pause; p.sol.ds=0.05; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,6);
%% branch plot
f=3; c=[1 9]; mclf(f); % c=[1 8]; % for V, which is not yet correct! 
plotbra('PH','pt20',f,c,'cl','k','lab',16,'bplab',1); 
plotbra('PHb','pt16',f,c,'cl',p2pc('gr1'),'bplab',1,'lab',15);
plotbra('za','pt6',f,c,'cl',p2pc('g1'),'lab',6); 
plotbra('zb','pt9',f,c,'cl',p2pc('g2'), 'lab',6); 
plotbra('za2','pt6',f,c,'cl',p2pc('o1'),'lab',6); 
plotbra('zb2','pt7',f,c,'cl',p2pc('o2'), 'lab',60); 
ylabel('A');
%% soln plots 
p2pglob.vi=[20,30]; p2pglob.showbd=2; p2pglob.edc='none'; p2pglob.cb=1; 
%%
plotsol('PH','bpt1'); pause; plotsol('PH','pt16'); pause 
plotsol('PHb','bpt1'); pause; plotsol('PHb','pt15');
%%
plotsol('za','pt6'); pause; plotsol('zb','pt6'); pause; 
plotsol('za2','pt6'); 
%%
ploteigsys('PH','bpt1',3); 
%%
p=loadp('PHb1','pt8','PHb1r'); plotsol(p); p2pglob.pbctol=1e-3; 
p.sw.rlong=1; sig=0.1; p.sw.nobdref=0; p=refineX(p,sig); 
p.sol.ds=0.01; p.nc.tol=1e-3; 
%% qswibra; ABE work, but looks like just the 2 orientations come back 
aux.besw=1; aux.m=2; aux.ral=1; p1=qswibra('PH','bpt1',aux);
%% note: ds=-0.05 jut flips axis! 
p=seltau(p1,2,'mm',2); p.sol.ds=0.05; p.sw.bifcheck=0; p=cont(p,10);
%% continuation alternating with moveX, and refine and degcoarsen, parameters:
p=loadp('PH1r','pt22','du');  ds=0.01; p.sol.ds=ds; p.nc.dsmax=5*ds; 
nis=5; ncs=1; dt=0.05; nit=1; sigr=0.01; sigc=0.01;  cit=5; keepbd=1; p.nc.tol=1e-3; 
p2pglob.pbctol=1e-2; 
%%
for i=1:4;     % outer loop, 
  for j=1:nis;  % inner loop, alternate moveX and cont  
      p=moveX(p,dt,nit); pplot(p,20); p=cont(p,ncs); 
  end  
p=refineX(p,sigr); % 
p=degcoarsenX(p,sigc,cit,keepbd);  
end


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
p2pglob.pbctol=1e-10; p=loadp('Pb','pt15','Pbc'); pplot(p,1); nit=2; keepbd=1; 
sig=0.1; p=degcoarsenX(p,sig,nit,keepbd); p.sol.ds=-p.sol.ds; pause; p=cont(p,5); 
%% repeat degcoarsen, OK 
p2pglob.pbctol=1e-10; p=loadp('Pbc','pt20','Pbc'); sig=0.05; 
p=degcoarsenX(p,sig,nit,keepbd); pause; p=cont(p,10); 