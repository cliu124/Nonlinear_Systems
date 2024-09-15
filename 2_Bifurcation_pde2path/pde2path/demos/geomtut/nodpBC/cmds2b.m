%% double cell, 2nd BP,  trying to solve quadratic BE
aux.besw=1; aux.ali=[1 3];  aux.isotol=1e-4; 
aux.mu2=0.2; aux.m=4; aux.ral=1; p1=qswibra('lN','bpt2',aux);
%% same as lN2a from cmds2 
p=seltau(p1,1,'q1',2); p.nc.tol=1e-4; p.sol.ds=-0.01;p.sw.bifcheck=0; p=cont(p,2); 
%% switch on PC and cont further 
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=-0.01; p=cont(p,15); 
%% again same as lN2a from cmds2 
p=seltau(p1,2,'q2',2); p.nc.tol=1e-4; p.sol.ds=-0.01;p.sw.bifcheck=0; p=cont(p,2); 
%% switch on PC and cont further 
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=-0.01; p=cont(p,15); 
%% --- 2nd BP, trying to solve cubic BE; finds no non-isol. solns at all 
aux.besw=1; aux.ali=[1 3];  aux.isotol=1e-4; 
aux.mu2=0.2; aux.m=4; aux.ral=1; p1=cswibra('lN','bpt2',aux);
%% ----------------------- other h0; -mainly so see how BP2 from h0=1 separates 
close all; keep pphome; clear classes; 
%% 
global p2pglob; p2pglob.vi=[10 20]; p2pglob.edc='k'; p2pglob.faceal=1; p2pglob.cut=1;  
p2pglob.showN=0; p2pglob.tsw=3; p2pglob.pbctol=1e-8; 
h0=1; dpre='./';  % output, uncomment as desired; first org 
h0=1.1; dpre='p2/';
h0=0.8; dpre='p1/'; 
%%
V0=0; A0=0; a=2.4; r=0.2; del=1; sx=0; sy=0; sz=0; sphi=0;
par=[h0; V0; A0; a; r; del; sx; sy; sz; sphi]; % del=height, 10=rot
%     1             5       7           10                                                   
lx=2*pi; ly=pi; sym=0; p=[]; ny=70; nx=60; %start with rather coarse mesh
p=nodbuckinit(p,lx,ly,nx,ny,par,sym); p=setfn(p,[dpre 'lN']); huclean(p);
pplot(p); p.sw.Xcont=1; p0=p;
%% 1st step  
p=p0; p.nc.ilam=[6 7 8 9]; p.nc.nq=3; p.nc.tol=1e-6; mclf(2);  
p.fuha.qf=@qf; p.fuha.qfder=@qjac; p.sw.qjac=1; p.sw.jac=0;
p.sol.ds=-0.01; p.nc.dsmax=0.04; p.sw.para=2; p.plot.bpcmp=19; p.nc.neig=50; 
p.nc.mu1=0.2; p.nc.mu2=0.05; p.nc.foldtol=1e-2; p.nc.bisecmax=6; p=cont(p,1);
%% go 
p.nc.tol=1e-8;p.sol.ds=-0.01;p.nc.dsmax=0.05; p.nc.dsmin=1e-4;
p=cont(p,30); 
%% other direction 
p=loadp([dpre 'lN'],'pt1',[dpre 'lNb']); p.sol.ds=-p.sol.ds; p=cont(p,8); 
%% step through BPs and check evals and kernels; 
% h0=1; dubious BP is bpt2;  h0=1.1; dubious BP is bpt2; h0=0.8; dubious BP is bpt1
% but Evals clearly better separated at h0=0.8, h0=1.1 than at  h0=1; 
aux.besw=0; aux.mu2=0.2; aux.m=4; aux.ral=1; p1=qswibra([dpre 'lN'],'bpt2',aux);
%% remainder not really relevant. 
p=gentau(p1,1,'p2/lN1a',2); p.sol.ds=0.1; p.nc.tol=1e-4; p.sw.bifcheck=0; 
p=cont(p,2); % two steps without rot.PC; then switch on and cont further
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=-0.01; p=cont(p,15); 
%% branch plot
f=3; mclf(f); p2pglob.pbctol=1e-4;  ylab='r'; xlab='\delta'; 
xlab='\delta';  c=[6 19]; 
plotbra([dpre 'lN'],'pt30',f,c,'cl','k'); 
xlabel(xlab); ylabel(ylab);  
%% soln plot
p2pglob.vi=[30,50]; p2pglob.tsw=6; p2pglob.edc='k';  
p2pglob.cut=1; plotsol([dpre 'lN'],'bpt1'); 
