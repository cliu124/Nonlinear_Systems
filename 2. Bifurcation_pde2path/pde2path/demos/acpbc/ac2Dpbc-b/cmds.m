%% ac2Dpbc cmds, square, pBC in x and y 
close all; keep pphome; 
%% init, parameters [c lambda gamma sx=x-speed sy=y-speed] 
p=[]; par=[1 -0.3 1 0 0]; lx=pi; ly=pi; p=acinit(p,lx,ly,20,par); 
p.nc.lammax=4; p.sol.ds=0.1; p.nc.dsmax=0.1; p=setfn(p,'tr'); p.sw.jac=1; 
%% first continue in d
p.nc.ilam=2; p=cont(p,30); 
%% 1st BP is simple
p=swibra('tr','bpt1','hom', 0.1); p=cont(p,40);
%% 3rd BP is 'effectively' simple, but needs 1 PC 
p=swibra('tr','bpt3','3s', 0.1); p=cont(p,2); p=qxon(p); p=cont(p,30); 
%% 2nd BP, stripes and spots, first only compute and plot kernel
aux=[]; aux.besw=0; p0=cswibra('tr','bpt2',aux); 
%% Now select 'active kernel list' aux.ali, derive and solve CBEs 
aux.ali=[1 2 3];aux.soltol=5e-10;aux.hasker=1;aux.besw=1; p0=cswibra(p0,aux); 
%% horizontal stripes, 3 initial steps, then switch on PC in y 
p=seltau(p0,1,'2hs',3); p.sol.ds=-0.05; p=cont(p,2); p=qyon(p); p=cont(p,30); 
%% vertical stripes 
p=seltau(p0,2,'2vs',3); p.sol.ds=-0.05; p=cont(p,2); p=qxon(p); p=cont(p,30); 
%% spots, need both PCs 
p=seltau(p0,8,'2sp',3); p.sol.ds=-0.1; p=cont(p,1); p=qxyon(p); p=cont(p,30); 
%% plot BD 
f=3; c=0; figure(f); clf;
plotbra('tr',f,c,'cl','k','bplab',[2 3],'lsw',0); 
plotbra('hom',f,c,'cl','b','lsw',0); 
plotbra('2hs',f,c,'cl','r','lab',14); 
plotbra('2sp',f,c,'cl','m','lab',13); 
plotbra('3s',f,c,'cl','b','lsw',0); 
axis([-0.3 3 0 7]); 
%% solution plots 
plotsol('2hs','pt14'); noticks; pause; plotsol('2sp','pt13'); noticks;
%% continuation in gamma of FPT1 on 2sp: (for checking fold-cont with contraints) 
p=[]; figure(2); clf; p=spcontini('2sp','fpt1',3,'2spfc'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; 
p.nc.del=1e-4; p.sw.spjac=1; p.fuha.spjac=@spjac; 
p.sw.spqjac=1; p.fuha.spqjac=@stanspqjac; 
p.nc.tol=1e-6; huclean(p); %[Ja, Jn]=spjaccheck(p); pause 
tic; p=cont(p,10); toc
%% switch back to regular continuation from one of the fold points
p=spcontexit('2spfc','pt10','2sp2-a'); p=resetc(p); p.nc.dsmax=0.1; p.sw.bifcheck=0; 
p.plot.bpcmp=0; p.nc.lammin=-5; p.sol.ds=1e-3; clf(2); huclean(p); p=cont(p,10); 
p=loadp('2sp2-a','pt5','2sp2-b'); p.sol.ds=-p.sol.ds; p=cont(p,20); % other direction 
%% plot BD of the two 2sp branches 
f=3; c=0; figure(f); clf;
plotbra('2sp',f,c,'cl','m','lsw',0); 
plotbra('2sp2-a','pt15',f,c,'cl',p2pc('r1'),'lsw',0); 
plotbra('2sp2-b',f,c,'cl',p2pc('r2'),'lsw',0,'fp',6); 
axis([0.7 1.5 0 6]); 