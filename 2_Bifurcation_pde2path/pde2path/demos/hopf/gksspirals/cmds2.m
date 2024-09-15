%% script for bif and cont of spirals ([Hotheo, §3.2.2]. RW spirals computed 
% (AND CONTINUED IN DELTA (dom-size)) as RELATIVE EQUILIBRIA. for old version see cmds2_old 
% See also cmds3.m, which also deals with mRWs, i.e., meandering spirals
close all; format compact; keep pphome;  
%% init, and find bifpoints from trivial branch, same HBPs as in cmds1
% but with al=1, thus recompute bifpoints to have the correct NF coeff
ndim=2; dir='tr-b'; p=[]; lx=pi; nx=35; r=-0.25; al=1; del=1; s=0; par=[r al del s]; 
p=rotinit(p,nx,par); p=setfn(p,dir); p.nc.neig=50; 
p.sol.ds=0.01; p.nc.dsmax=0.01; p.nc.ilam=1; p.nc.bisecmax=20; p=cont(p,60); 
%% RWs via twswibra 
para=4; ds=0.4; dsmax=1.7; xi=1e-2; nsteps=50; lammax=3.1; 
figure(2); clf; aux=[]; usrlam=[0 1 2 3]; spar=4; 
for bnr=2; %[2 3 5 6 7] 
switch bnr
 case 2; kwnr=1; p=twswibra('tr-b','hpt2',spar, kwnr,'rs2',aux); 
 case 3; kwnr=2; p=twswibra('tr-b','hpt3',spar, kwnr,'rs3',aux); 
 case 5; kwnr=3; p=twswibra('tr-b','hpt5',spar, kwnr,'rs5',aux); 
 case 6; kwnr=1; p=twswibra('tr-b','hpt6',spar, kwnr,'rs6',aux); 
 case 7; kwnr=4; p=twswibra('tr-b','hpt7',spar, kwnr,'rs7',aux); 
end 
p.nc.dsmax=dsmax; p.nc.lammax=lammax; 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; 
p.u0x=p.mat.M\(p.mat.Krot*p.u0); p.u(1:p.nu)=p.u(1:p.nu)+0.1*p.tau(1:p.nu); 
plotsolu(p,p.u0,6,1,3); plotsolu(p,5*p.u0x,11,1,3); p.sw.bifcheck=2; 
p.nc.nq=1; p.nc.ilam=[1;spar];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac;
p.usrlam=usrlam; 
p.nc.tol=1e-4; p=cont(p,2); p.nc.tol=1e-8; p=cont(p,20); 
end 
%% continue in diff constant (domain size)
p=swiparf('rs2','pt14','rs2d',[3 4]); clf(2); p.sol.ds=-0.05; p.nc.dsmax=0.21; 
p.usrlam=[0.2 0.3 0.4 0.5]; p.nc.lammin=0.1; p.nc.neig=30; p=cont(p,40);
%% BD of cont in delta
bpcmp=8; figure(3); clf; plotbra('rs2d','pt20',3,bpcmp,'cl','r','lab',[14,19]); 
xlabel('\delta'); ylabel('||u||_*'); 
%% continue in r again
p=swiparf('rs2d','pt19','rs2dr',[1 4]); clf(2); p.sw.verb=2; p.sol.ds=-0.05; 
p.usrlam=[0 1 2 3]; p.nc.lammin=-0.25; p=cont(p,40);
%% BD of cont in r
bpcmp=8; figure(3); clf; plotbraf('rs2dr','pt30',3,bpcmp,'lab',[20 25],'cl','b'); 
xlabel('r'); ylabel('||u||_*'); 
%% BD of six branches, L^2
bpcmp=8; figure(3); clf; plotbra('tr',3,bpcmp,'cl','k'); 
plotbra('rs2','pt14',3,bpcmp,'lab',[10,14],'cl','r'); 
plotbra('rs3',3,bpcmp,'cl','m'); plotbra('rs5',3,bpcmp,'cl','b'); 
plotbra('rs6',3,bpcmp,'cl','r'); plotbra('rs7',3,bpcmp,'cl','m'); 
axis([-0.3 3 0 1.2]); xlabel('r'); ylabel('||u||_*'); 
%% profiles from cont in r 
plotsol('rs2','pt10',1,1,2); clplot(1); pause; plotsol('rs2','pt14',1,1,2); clplot(1); pause 
plotsol('rs3','pt15',1,1,2); clplot(1); pause; plotsol('rs5','pt16',1,1,2); clplot(1); pause 
plotsol('rs6','pt12',1,1,2); clplot(1); pause; plotsol('rs7','pt17',1,1,2); clplot(1); 
%% profiles from cont in delta
plotsol('rs2d','pt14',1,1,2); clplot(1); pause; plotsol('rs2d','pt19',1,1,2); clplot(1); 
%% profiles for subsequent cont in r
plotsol('rs2dr','pt20',1,1,2); clplot(1); 
%% spectral plot: (at point where rs2 gains stability) 
p=loadp('rs2','hpt1'); wnr=10; p.nc.neig=100; [muv,ev]=plotspec(p,wnr); 