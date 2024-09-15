%% two inner buckles nodoid, i.e., twice the minimal cell of cmds1
close all; keep pphome; clear classes; % modification in stanpdeo2D 
%% 
global p2pglob; p2pglob.vi=[10 20]; p2pglob.edc='k'; p2pglob.faceal=1; p2pglob.cut=1;  
p2pglob.showN=0; p2pglob.tsw=3; p2pglob.pbctol=1e-8; 
V0=0; A0=0; h0=1; a=2.4; r=0.2; del=1; sx=0; sy=0; sz=0; sphi=0;
par=[h0; V0; A0; a; r; del; sx; sy; sz; sphi]; % del=height, 10=rot
%     1             5       7           10                                                   
lx=2*pi; ly=pi; sym=0; p=[]; ny=70; nx=60; %start with rather coarse mesh
p=nodbuckinit(p,lx,ly,nx,ny,par,sym); p=setfn(p,'lN'); huclean(p);
pplot(p); p.sw.Xcont=1; p0=p;
%% 1st step  
p=p0; p.nc.ilam=[6 7 8 9]; p.nc.nq=3; p.nc.tol=1e-6; mclf(2);  
p.fuha.qf=@qf; p.fuha.qfder=@qjac; p.sw.qjac=1; p.sw.jac=0;
p.sol.ds=-0.01; p.nc.dsmax=0.04; p.sw.para=2; p.plot.bpcmp=19;
p.nc.mu1=0.2; p.nc.mu2=0.05; p.nc.foldtol=1e-2; p.nc.bisecmax=10; p=cont(p,1);
%% refine initial mesh 
p2pglob.pbctol=1e-3; sig=0.5; p=loadpp('lN','pt1'); p.sw.rlong=1; p.sw.nobdref=0; 
p=refineX(p,sig); p2pglob.pbctol=1e-3; %sig=0.25; p=refineX(p,sig); 
%% go 
p.nc.tol=1e-8;p.sol.ds=-0.01;p.nc.dsmax=0.05; p.nc.dsmin=1e-4;
p=cont(p,30); 
%%
p=swibra('lN','bpt1','lN1');p.sol.ds=-0.1; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,2);
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1;p.nc.tol=1e-4; p.sol.ds=-0.01; p=cont(p,15);
%% 2nd BP possibly double, but algebraic bifurcation equations (ABEs) only yield 
% pure modes (see cmds2b). Hence, here switch off ABEs for speed. 
% For small ds, two separate BPs may also be detected; but treating them
% via the qswibra-trick should always work 
aux.besw=0; % switch off ABE computations, just compute kernel
aux.mu2=0.01; aux.m=4; aux.ali=[]; p1=qswibra('lN','bpt2',aux);
%% go on 1st 
p=gentau(p1,1,'lN2a',2); p.sol.ds=0.1; p.nc.tol=1e-4; p.sw.bifcheck=0; 
p=cont(p,2); % two steps without rotational PC; then switch on and cont further
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-4; p.sol.ds=-0.01; p=cont(p,15); 
%% go on 2nd 
p=gentau(p1,[0 0 1],'lN2b',2); p.sol.ds=0.025; p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4);
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=-0.01; p=cont(p,15); 
%% depending on ds, the BP at del\approx 1.1 may be BP3 or BP4; adapt accordingly 
p=swibra('lN','bpt4','lN3');p.sol.ds=-0.1; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,2);
p.nc.nq=4; p.nc.ilam=[6 7 8 9 10]; p.fuha.qf=@qfrot; p.fuha.qfder=@qjacrot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=-0.01; p=cont(p,15); 
%% branch plot
f=3; mclf(f); p2pglob.pbctol=1e-4;  ylab='r'; xlab='\delta';  c=[6 19]; 
plotbra('lN','pt30',f,c,'cl','k','labi',[0],'lp',30);
plotbra('lN1','pt20',f,c,'cl','b','lab',[15]); 
plotbra('lN2a',f,c,'cl','r','lab',[10 18]); 
plotbra('lN2b','pt22',f,c,'cl',p2pc('g1'),'lab',[10 22],'fp',1); 
plotbra('lN3','pt8',f,c,'cl','m','lab',[8]); 
xlabel(xlab); ylabel(ylab);  
%% soln plot
p2pglob.vi=[30,50]; p2pglob.tsw=6; p2pglob.edc='k';  p2pglob.cb=0; p2pglob.showN=1; 
p2pglob.cut=1; plotsol('lN','bpt1'); pause; 
p2pglob.cut=1; plotsol('lN1','pt15'); pause;
p2pglob.showbd=2; 
%%
plotsol('lN2a','pt10');pause; plotsol('lN2a','pt18');pause;p2pglob.cut=0;  plotsol('lN2a','pt18'); pause; 
p2pglob.cut=1; plotsol('lN2b','pt10'); pause; plotsol('lN2b','pt22'); pause; 
p2pglob.cut=0; plotsol('lN2b','pt22'); pause; plotsol('lN3','pt8'); 