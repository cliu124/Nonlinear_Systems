%% Ac on equil.triange; D_3 symmetry 
close all; keep pphome; 
%% init, and localization of 2n BP
p=[]; par=[0.1 1.5 0.2 0]; % parameters [c lambda cub al]
hmax=0.07; r=1; p=acinit(p,r,hmax,par); p.np, pause 
p.nc.ilam=2; p.nc.lammax=1.9; p.sol.ds=0.02; p.nc.dsmax=0.02; p=setfn(p,'tr'); 
p.usrlam=[1.88];  p=findbif(p,1);
%% 2nd BP, double, Z2 branches alway (any al) bifurcate 
% Z3 branches only exists for al=0 (pitchforks) 
aux=[]; aux.m=2; p0=qswibra('tr','bpt1',aux); 
%% Z2 branches, backward and forward
p=gentau(p0,1); p.sol.ds=0.01; p=setfn(p,'b2ta'); p=cont(p,40);
p=gentau(p0,1); p.sol.ds=-0.01; p=setfn(p,'b2tab'); p=cont(p,40);
%% Z3 branches, only for al=0, and then pitchforks 
p=gentau(p0,[0 1]); p.sol.ds=0.03; p=setfn(p,'b2pa'); p=cont(p,40);
%% BD for al=0, ordinates are corner values 
f=3; c=9; figure(f); clf;
plotbra('b2ta','pt32',f,9,'cl','b','lsw',0,'lab',29); 
plotbra('b2tab','pt32',f,9,'cl','b','lsw',0);
plotbra('b2tab','pt32',f,8,'cl',p2pc('b1'),'lsw',0,'lab',29);
plotbra('b2ta','pt32',f,8,'cl',p2pc('b1'),'lsw',0); 
plotbra('b2pa','pt28',f,7,'cl',p2pc('r1'),'lsw',0); 
plotbra('b2pa','pt28',f,8,'cl',p2pc('r1'),'lsw',0); 
plotbra('b2pa','pt28',f,9,'cl',p2pc('r3'),'lsw',0,'lab',25); 
axis([1.75 1.9 -1.3 1.3]); xlabel('\lambda'); ylabel('u(x_0)'); box on; 
%% soln-plots (al=0) 
v=[-40,40]; v=[0,90];
plotsol('b2ta','pt29',1,1,2);view(v); pause
plotsol('b2pa','pt25',1,1,2);view(v);
%% Now al\ne 0; first continue phi2-branch at lam approx 0.185 in al
p=swiparf('b2pa','pt20','b2paq1',4); p.sol.ds=0.001; p=cont(p,20); 
%% continue back in lam again (al approx 1.5e-2) to find imperfect bifs 
p=swiparf('b2paq1','pt5','b2pa2',2); p.sol.ds=-0.001; p=cont(p,30); 
%% other direction 
p=swiparf('b2paq1','pt5','b2pa2b',2); p.sol.ds=0.001; p=cont(p,30); 
%% continue phi1 branch at al\ne 0  (choose al=0.2) for graphical reasons 
p=gentau(p0,1); p.sol.ds=0.01; p.u(p.nu+4)=0.2; p=setfn(p,'b2ta2'); p=cont(p,40);
%% other direction 
p=gentau(p0,1); p.sol.ds=-0.01; p.u(p.nu+4)=0.2; p=setfn(p,'b2ta2b'); p=cont(p,40);
%% BD for al\ne 0
f=4; c=9; figure(f); clf;
plotbra('b2ta2','pt41',f,8,'cl','b','lsw',0); 
plotbra('b2ta2b','pt41',f,8,'cl','b','lsw',0);
plotbra('b2ta2','pt41',f,9,'cl',p2pc('b3'),'lsw',0,'lab',38); 
plotbra('b2ta2b','pt40',f,9,'cl',p2pc('b3'),'lsw',0,'lab',20); 
plotbra('b2pa2','pt23',f,7,'cl',p2pc('r1'),'lsw',0); 
plotbra('b2pa2b','pt14',f,7,'cl',p2pc('r1'),'lsw',0);
plotbra('b2pa2','pt23',f,8,'cl',p2pc('r2'),'lsw',0);
plotbra('b2pa2b','pt14',f,8,'cl',p2pc('r2'),'lsw',0);
plotbra('b2pa2','pt23',f,9,'cl',p2pc('r3'));
plotbra('b2pa2b','pt14',f,9,'cl',p2pc('r3'));
axis([1.74 1.9 -1 1.8]); xlabel('\lambda'); ylabel('u(x_0)'); box on; 
%% soln plots 
plotsol('b2ta2','pt38',1,1,2);view(v);pause 
plotsol('b2ta2b','pt20',1,1,2);view(v);pause 
plotsol('b2pa2','pt19',1,1,2);view(v);pause 
plotsol('b2pa2b','pt10',1,1,2);view(v);
%% secondary bifs 
p=swibra('b2ta','bpt1','b2a1',0.1); p.pm.resfac=1e-6; p=pmcont(p,20); 
%%
p=swibra('b2pa','bpt1','b2p1',0.02); p.pm.resfac=1e-6; p=pmcont(p,20);