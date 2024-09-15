%% One inner buckle nodoid
close all; keep pphome; clear classes; % modification in stanpdeo2D 
%% 
global p2pglob; p2pglob.vi=[10 40]; p2pglob.edc='k'; p2pglob.faceal=1; p2pglob.cut=0;  
p2pglob.tsw=3; p2pglob.pbctol=1e-3; % weak tol to identify boundaries 
V0=0; A0=0; h0=1; a=2.4; r=0.2; del=1; sx=0; sy=0; sz=0; sphi=0;
par=[h0; V0; A0; a; r; del; sx; sy; sz; sphi]; % del=height, 10=rot
lx=pi; ly=pi; sym=0; p=[]; ny=70; nx=30; %start with rather coarse mesh
p=nodbuckinit(p,lx,ly,nx,ny,par,sym); p=setfn(p,'bN'); pplot(p); p.sw.Xcont=1; 
%% 1st step  
p.nc.ilam=[6 7 8 9]; p.nc.nq=3; % 3 translational constraints 
p.fuha.qf=@qf; p.fuha.qfder=@qjac; p.sw.qjac=1; p.sol.ds=-0.01; p=cont(p,1);
%% refine initial mesh (twice), in particular at boundary 
sig=0.2; p=loadpp('bN','pt1'); p.sw.rlong=1; p.sw.nobdref=0; 
p=refineX(p,sig); sig=0.25; p=refineX(p,sig); p=retrigX(p); 
%% go further 
p.nc.tol=1e-8;p.sol.ds=-0.01;p.nc.dsmax=0.04; p.nc.dsmin=1e-6; p=cont(p,30); 
%%
p=swibra('bN','bpt1','bN1');p.sol.ds=-0.1; p.nc.tol=1e-6; pause; p.sw.bifcheck=0; p=cont(p,2);
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1;p.nc.tol=1e-4; p.sol.ds=-0.01; p=cont(p,15);
%%
p2pglob.pbctol=1e-3; aux.besw=0; aux.mu2=0.2; aux.m=2; aux.ral=1; p1=qswibra('bN','bpt2',aux);
p=gentau(p1,[0.5 0.5],'bN2',2); p.sol.ds=0.01; p.nc.tol=1e-3; 
p.nc.dsmax=0.1; p.sw.bifcheck=0; p=cont(p,2);
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=-0.01; p=cont(p,15); 
%% other direction 
p2pglob.pbctol=1e-2; p=loadp('bN','pt1','bNb'); p.sol.ds=-p.sol.ds; p.nc.tol=1e-6; p=cont(p,10);
 %% branch plot
f=3; mclf(f); p2pglob.pbctol=1e-4;  ylab='r'; xlab='\delta'; 
xlab='\delta';  c=[6 19]; %xlab='z'; ylab='A'; 
plotbra('bN','pt30',f,c,'cl','k','lab',20);
plotbra('bNb','pt10',f,c,'cl',p2pc('gr1'),'lab',10);
plotbra('bN1','pt22',f,c,'cl','b','lab',[15 22]); 
plotbra('bN2','pt20',f,c,'cl','r','lab',[20]); 
%plotbra('bN1r','pt20',f,c,'cl','b','lab',[20]); 
%plotbra('bN2r','pt20',f,c,'cl','r','lab',[20]); 
xlabel(xlab); ylabel(ylab);  
%% soln plot
p2pglob.vi=[30,50]; p2pglob.tsw=5; p2pglob.edc='k'; p2pglob.cut=0; 
plotsol('bNb','pt10'); pause; plotsol('bN','pt3'); pause; plotsol('bN','pt20'); pause 
plotsol('bN1','pt15');pause;  plotsol('bN1','p20');pause; plotsol('bN2','pt15'); pause 
plotsol('bN1r','pt22'); pause; plotsol('bN2r','pt20');







