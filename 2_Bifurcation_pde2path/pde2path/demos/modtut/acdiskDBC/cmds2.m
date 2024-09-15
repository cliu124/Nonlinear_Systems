%% ac2D disk, Nonhom. DBCs, 
close all; keep pphome; 
%% init, and cont of prim branch 
p=[]; s=0; par=[0.5 -0.1 0 -1 s 0.1]; % c,lam,quad,cubic; coeff of bdry term 
lx=5; ly=1; nr=16; nphi=32; p=acinit(p,lx,nr,nphi,par); p=setfn(p,'i1'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=0.5; p.sw.verb=2; p.nc.neig=20; p.nc.lammax=1; p.nc.eigref=-1;  
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8; p.plot.bpcmp=7;  p=cont(p,10);
%% imperfect bif to i2 
p=swibra('i1','bpt1','i2'); p.sw.verb=2; pause; p.sw.bifcheck=0; p=cont(p,10); 
%% init at larger lam to get other branch 
p=[]; s=0; par=[0.5 0.5 0 -1 s 0.1]; % c,lam,quad,cubic; coeff of bdry term 
lx=5; ly=1; nr=16; nphi=32; p=acinit(p,lx,nr,nphi,par); p=setfn(p,'i3'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=0.5; p.sw.verb=2; p.nc.neig=20; p.nc.lammax=1; p.nc.eigref=-1;  
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8; p.plot.bpcmp=7;  p=cont(p,10);
%% other direction 
p=loadp('i3','pt0','i3b'); p.sol.ds=-0.1; p=cont(p,20); 
%% cont i2 in bd-amp 
p=swiparf('i2','pt6','i2ga',6); p.nc.dsmax=0.1; p.nc.lammax=2; p=cont(p,70); 
%% BD plot (cont in lam) 
f=3; c=7; figure(f); clf;
plotbra('i1','pt10',f,c,'cl','k','lab',9); 
plotbra('i2',f,c,'cl','b','lab',6); 
plotbra('i3','pt10',f,c,'cl','r','lab',9); plotbra('i3b',f,c,'cl','m','lab',12); 
xlabel('\lambda'); ylabel('||u||_2');
axis([-0.1 1.01 0 0.9]); 
%% BD plot (cont in ga) 
f=3; c=7; figure(f); clf;
plotbra('i2ga','pt70',f,c,'cl','b','lab',[0, 20 70]); 
xlabel('\gamma'); ylabel('||u||_2');
%% soln plots 
plotsol('i1','pt9'); pause; plotsol('i2','pt6'); pause; plotsol('i2ga','pt20'); 
pause; plotsol('i2ga','pt70'); pause; plotsol('i3','pt9'); pause; plotsol('i3b','pt12'); 