%% cont in lam1, at c0=0.7; 4 non-axi branches bifurcating at lam1<0, but 
% continuing to lam1>0 
close all; keep pphome; 
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.cb=1;
%% start with Helfrich cylinder but different mesh compared to cmds1, cont in c0
al=1; lx=1; ly=pi; nx=35; ny=round(2*al*nx);  sym=0; p=[]; l1=0.25; c0=0; srot=0; 
par=[al; l1; c0; srot]; p=cylinit(p,lx,nx,ny,par,sym); p=setfn(p,'l00'); pplot(p,1); 
p.nc.dsmax=0.2; p.nc.eigref=-50; p.nc.mu1=5; p.nc.mu2=0.2; p.nc.bisecmax=12; p.np
pause; p=cont(p,10); 
%% cont in lam1=tension, c0=0.7 
mclf(2); p=swiparf('l00','pt5','l0',2); p.nc.dsmax=0.2; p.usrlam=[]; 
p.sw.cdbb=1; p.sol.ds=-0.1; p.nc.tol=1e-6; tic; p=cont(p,20); toc
%% other direction 
p=loadp('l0','pt0','l0b'); p.sol.ds=-p.sol.ds; p=cont(p,10); 
%% 1st non-axisym branch 
p=swibra('l0','bpt1','l1',0.1); pause; p.nc.dsmax=0.2; p.nc.tol=1e-4; 
p.sw.bifcheck=0; p=cont(p,5); 
%% 1st non-axisym branch via qswibra 
aux=[]; aux.m=2; aux.besw=0; p0=qswibra('l0','bpt1',aux); pause;  
p=gentau(p0,[0 1],'l1',2); p.sol.ds=0.1; 
p.nc.dsmax=0.2; p.nc.tol=1e-3; 
p.sw.bifcheck=0; p=cont(p,5); 
%% switch on rotational PC, 
p=loadp('l1','pt4','l1q'); p.nc.nq=1; p.nc.ilam=[2 4]; 
p.fuha.qf=@qfrot; p.fuha.qfder=@qrotder; p.sw.qjac=1; % now in cylinit
p.nc.tol=1e-5; p=cont(p,5); getaux(p)' 
%% alternate ref. and cont
p=loadp('l1q','pt8','l1qr'); plotsol(p); sig=0.05; p.sol.ds=0.05; p.nc.dsmax=0.1; 
p.sw.spcalc=0; p.sw.bifcheck=0; 
for i=1:5; p=refineX(p,sig); p=cont(p,4); end; getaux(p)' 
%% 2nd non-axisym branch 
p=swibra('l0','bpt3','l2',0.1); pause; p.nc.dsmax=0.2; p.nc.tol=1e-4; 
p.sw.bifcheck=0; p=cont(p,5); 
%% switch on rotational PC, 
p=loadp('l2','pt3','l2q'); p.nc.nq=1; p.nc.ilam=[2 4]; p.nc.tol=1e-6; p=cont(p,6); getaux(p)' 
%% alternate ref. and cont
p=loadp('l2q','pt8','l2qr'); plotsol(p); sig=0.05;  p.sol.ds=0.05; p.nc.dsmax=0.05; 
for i=1:5; p=refineX(p,sig); p=cont(p,6); end; getaux(p)' 
%% 3rd non-axisym branch 
p=swibra('l0','bpt5','l3',0.1); pause; p.nc.dsmax=0.2; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,5); 
%% switch on rotational PC, 
p=loadp('l3','pt4','l3q'); p.nc.nq=1; p.nc.ilam=[2 4]; p=cont(p,4); getaux(p)' 
%% alternate ref. and cont
p=loadp('l3q','pt8','l3qr'); plotsol(p); sig=0.075; p.sol.ds=0.05; p.nc.dsmax=0.05; 
for i=1:5; p=refineX(p,sig); p=cont(p,5); end; getaux(p)'
%% 4th non-axisym branch 
p=swibra('l0','bpt6','l4',0.1); pause; p.nc.dsmax=0.2; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,5); 
%% switch on rotational PC, 
p=loadp('l4','pt4','l4q'); p.nc.nq=1; p.nc.ilam=[2 4]; p=cont(p,3); getaux(p)' 
%% alternate ref. and cont
p=loadp('l4q','pt7','l4qr'); plotsol(p); sig=0.075; p.sol.ds=0.05; p.nc.dsmax=0.05; 
for i=1:5; p=refineX(p,sig); p=cont(p,7); end 
%% branch plot, 
f=3; mclf(f); xlab='\lambda_1'; c=5; ylab='A'; ax=[-1.9 1 12 30]; 
%c=7; ylab='E'; ax=[-1.9 1 -25 25]; 
%c=8; ylab='\delta_{mesh}'; ax=[-1.9 1 -25 25]; 
plotbra('l0','pt38',f,c,'cl','b','lab',36); 
plotbra('l0b','pt10',f,c,'cl',p2pc('b2'),'lab',5); 
plotbra('l1qr','pt26',f,c,'cl','r','lab',[8 26]); 
plotbra('l2qr',f,c,'cl','m','lab',36); plotbra('l3qr',f,c,'cl',p2pc('o1'),'lab',30); 
plotbra('l4qr',f,c,'cl',p2pc('g1'),'lab',34); 
xlabel(xlab); ylabel(ylab); axis(ax); grid on
%% soln plots
p2pglob.vi=[60,20]; p2pglob.edc='none'; p2pglob.cm='parula'; p2pglob.tsw=0; p2pglob.showbd=2; 
pplot('l0b','pt5');  pause; p2pglob.edc='k'; pplot('l0','pt36');  pause; 
pplot('l1q','pt8'); pause; pplot('l1qr','pt26'); pause; p2pglob.edc='none';
pplot('l2qr','pt36'); pause; pplot('l3qr','pt30'); pause; pplot('l4qr','pt34');
%%
pplot('l1qr','pt28'); pause; pplot('l1qrc','pt27');
%%
plotspec('l4qr','pt20',6); 
%% degcoarsening-cont-loop, doesn't quite work 
p=loadp('l1qr','pt20','l1qrc'); plotsol(p); sigc=0.15; sigr=0.05; p.nc.dsmax=0.1; p.nc.ds=0.1; 
for i=1:4; p=degcoarsenX(p,sigc,4,1); p=cont(p,4); p=refineX(p,sigr); p=cont(p,2); end 
%% branch plot
f=3; c=10; mclf(f); 
plotbra('l1qr','pt28',f,c,'cl','r','lab',[10],'fp',20); 
plotbra('l1qrc','pt25',f,c,'cl',p2pc('r2'),'lab',[32,36],'fp',20); 
ylabel('max(h/r)'); 