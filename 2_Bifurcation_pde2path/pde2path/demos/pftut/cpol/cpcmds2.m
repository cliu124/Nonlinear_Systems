%% cell polarization, Cussedu-Edelstein-Keshet-..-Madzvamuse 2018, rough mesh, incl.other orientations 
close all; keep pphome; 
%% init and trivial branch, simple case with rough mesh, nu=3200
p=[]; e=0.125; k0=0.1; ga=1; R=1; bpstyle=4; m=11; s=0; mc=0; 
cut=[-2 0 -2]; par=[R k0 e ga s m mc]; % s=speed for x-PC, 
lx=pi; del=5e-2; ly=pi/2-del; 
nx=12; ny=8; ref=1; h0=0.125; odir='h1'; 
%nx=8; ny=7; ref=0; h0=0.2; odir='h2'; % just for testing 
p=cpinit(p,lx,ly,nx,ny,par,ref,bpstyle,cut,h0); 
p=setfn(p,odir); userplot(p,21); 
%% homogen. branch, continue with mass constraint 
p.nc.tol=1e-6; p.nc.ilam=[6 7]; p.nc.nq=1; p=cont(p,10);
%% patterned branch, N-S-orientation, hence no rot PC needed 
p=swibra(odir,'bpt1','b1',0.01);  p.nc.tol=1e-6; p.sw.bifcheck=0; p.sw.para=2; % run with poor tol  
p.nc.dsmax=0.4; p.nc.intol=-5e-3; p=cont(p,20);% relax instab. tolerance  
p.nc.dsmax=0.025; p.nc.intol=0; p=cont(p,15); % tighten instab.tol, reduce max ds 
%% BD's 
fnr=3; figure(3); clf; c=8; plotbra('h1','pt20',fnr,c); 
plotbra('b1','pt35',fnr,c,'cl','b','lab',[5 20]); ylabel('max(u)'); pause 
fnr=3; figure(3); clf; c=10; plotbra('h1','pt20',fnr,c); 
plotbra('b1','pt35',fnr,c,'cl','b','lab',[5 20]); ylabel('max(w)'); pause 
fnr=3; figure(3); clf; c=7; plotbra('h1','pt20',fnr,c); 
plotbra('b1','pt35',fnr,c,'cl','b'); 
plotbra('b1c','pt30',fnr,c,'cl','r','lab',[5 20]); ylabel('\delta'); 
%% soln plots
plotsol('b1','pt5'); pause; plotsol('b1','pt20'); pause; plotsol('b1c','pt5'); pause; plotsol('b1c','pt20');
%% check the threefold bif 
odir='h0'; 
aux=[]; aux.m=5; aux.besw=0; p0=cswibra(odir,'bpt1',aux); pause 
%aux.hasker=1; aux.ali=[1 2 3]; aux.besw=1; aux.isotol=1e-2; p0=cswibra(p0,aux);
p0.sw.bifcheck=0; 
%% cont other orientation with PC, 
p=seltau(p0,1,'b1b',3); p.usrlam=[]; p.nc.tol=1e-6; p.nc.dsmax=0.4; 
p.nc.intol=-5e-3; p=conpc(p,-0.05,2,20); 
%% cont other orientation with PC, 
p=gentau(p0,[0 1],'b1c'); p.usrlam=[]; p.nc.tol=1e-6; p.nc.dsmax=0.4; p.sw.para=2; p=cont(p,2), pause
p.nc.intol=-1e-2; p=conpc(p,-0.01,2,28); 
%%  testing and aux. stuff ... 
[Gu, Gn]=jaccheck(p); Gd=abs(Gu-Gn); em=max(max(Gd)); 
figure(10); clf; Ge=(Gd>em/5); spy(Ge)
%% check that points from surf map to points on surf of ball 
figure(10); clf; p.p2.grid.plotFaces([],'FaceAlpha',1); hold on; mapogrid(p,1);
%% patterned branch, N-S-orientation, hence no rot PC needed 
p=swibra('h2','bpt1','a1',0.01);  p.nc.tol=1e-6; p.sw.bifcheck=0; p.sw.para=2; % run with poor tol  
p.nc.dsmax=0.4; p.nc.intol=-5e-3; p=cont(p,25);% relax instab. tolerance  
p.nc.dsmax=0.025; p.nc.intol=0; p=cont(p,15); % tighten instab.tol, reduce max ds 
