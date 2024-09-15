close all; keep pphome; 
%% commands for Schnakenberg on BCC; too hard for octave at the moment 
p=[]; kc=sqrt(sqrt(2)-1); lx=sqrt(2)*pi/kc; ly=lx; lz=lx; par=[3.25, 0, 60]; 
nx=10; sw.sym=0; sw.ref=1; dir='BCC2'; % np=6006, still OK
p=schnakinit(p,[lx,ly,lz],nx,par,sw); p.pdeo.grid.plotFaces; p.np, pause 
%p=setilup(p,1e-4,300); p.fuha.lss=@lssAMG; % for large np 
p.nc.neig=10; p.sw.verb=2; p.nc.tol=1e-6; p.plot.pstyle=3; 
p.pm.resfac=1e-4; p.sol.ds=-0.1; p=setfn(p,dir); p=cont(p,3); 
%% qswibra; 3 tubes as kernel 
aux=[]; %aux.isotol=1e-16; 
aux.m=3; p0=qswibra('BCC2','bpt1',aux); 
p0.nc.dsmin=0.05; p0.pm.mst=4;  p0.sw.bifcheck=0; p0.pm.resfac=1e-6; p0.pm.mst=10; 
p0.file.smod=5; p0.nc.tol=1e-6; p0.nc.lammin=2.2; p0.sw.foldcheck=0; p0.sw.spcalc=0; 
%% select BCC and cont in both directions 
p=seltau(p0,4,'BCCa2',2); p.sol.ds=0.02;  p.plot.pstyle=3; p=pmcont(p,20); 
%p=seltau(p0,4,'BCCb2',2); p.sol.ds=-0.02; p=pmcont(p,5); 
%% use gentau for tubes 
p=gentau(p0,[0 0 1],'BCCt'); p=pmcont(p,20); 
%% plot BD 
fnr=3; cmp=4; figure(fnr); clf; %plotbra('BCC1','pt6',fnr,cmp,'cl','k','lp',3); 
plotbra('BCCa2','pt20',fnr,cmp,'cl',p2pc('r1'),'lab',20);
plotbra('BCCb2','pt5',fnr,cmp,'cl',p2pc('r2'),'lab',5); 
plotbra('BCCt','pt20',fnr,cmp,'cl','b','lab',20); 
xlabel('\lambda'), ylabel('max(|u|)'); box on; 
%% solns plots 
v=[-40,20]; plotsol('BCCa2','pt20'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
%plotsol('BCCa','pt25'); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('BCCb2','pt5',1,1,2); view(v); xlabel(''); ylabel(''); zlabel(''); pause; 
plotsol('BCCt','pt20',1,1,2); view(v); xlabel(''); ylabel(''); zlabel(''); 