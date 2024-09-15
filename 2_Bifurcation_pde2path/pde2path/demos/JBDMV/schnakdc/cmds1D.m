close all; keep pphome; 
%% DC-Schnakenberg; don't choose m too large (non-conv at small lam) 
p=[]; par=[3.2 1 100 1.5]; nx=300; lx=12; % [lam c d m] 
p=schnakinit(p,lx,nx,par); p=setfn(p,'tr'); 
p.zc=1; p.zc2=0; % zero-correction in Newton and in sG 
p=findbif(p,6); p.sw.bifcheck=0; p.nc.dsmax=1; p=cont(p,20); 
%% 
p=swibra('tr','bpt1','T1',0.1); p=cont(p,50);
p=swibra('tr','bpt2','T2',0.1); p=cont(p,50);
p=swibra('tr','bpt3','T3',0.1); p=cont(p,50);
%% steady from branch with DC 
p=swibra('T1','bpt1','T1-1',0.05); p.nc.dsmax=0.1; p.nc.tol=2e-8; p=cont(p,120);
%% Hopf
aux.dlam=0; aux.tl=40; p=hoswibra('T1','hpt1',0.1,4,'T1-h1',aux); p.file.smod=2; 
p.nc.tol=1e-6; p.hopf.flcheck=0; p.fuha.outfu=@hobra; p.usrlam=[]; 
p.x0i=p.np; p=setbel(p,1,1e-6,10,@lss); p=cont(p,10); 
%% Hopf
aux.dlam=0; aux.tl=50; p=hoswibra('T3','hpt1',0.1,4,'T3-h1',aux); p.file.smod=2; 
p.nc.tol=1e-6; p.hopf.flcheck=0; p.fuha.outfu=@hobra; p.usrlam=[]; 
p.x0i=p.np; p=setbel(p,1,1e-6,10,@lss); p=cont(p,10); 
%% plot BD 
fnr=3; figure(fnr); clf; c=6; plotbra('tr','pt50',fnr,c,'cl','k','fp',1); 
plotbra('T1',fnr,c,'cl','b', 'lab',25); 
plotbra('T1-H1','pt8',fnr,c,'cl',p2pc('o1'),'lab',[2 8]); 
plotbra('T2',fnr,c,'cl',p2pc('b1'),'lab',30); 
plotbra('T3',fnr,c,'cl',p2pc('g1'),'lab',40); 
plotbra('T3-2',fnr,c,'cl',p2pc('r2'),'lab',15); 
plotbra('T1-1','pt140',fnr,c,'cl',p2pc('r3'), 'lab',[40 100],'lp',105); 
plotbra('T3-H1','pt8',fnr,c,'cl',p2pc('o2'),'lab',[8]); 
ylabel('max(u)');
%% soln plot, steady 
mclf(1); 
mypsol('T1','bpt1'); pause; mypsol('T1','hpt1'); pause; mypsol('T1','pt25'); pause; 
mypsol('T1-1','pt40'); pause; mypsol('T1-1','pt100'); pause; 
mypsol('T3','bpt2'); pause; mypsol('T3','hpt1'); pause; mypsol('T3','pt40'); pause;
mypsol('T2','bpt1'); pause; mypsol('T2','hpt1'); pause; mypsol('T2','pt30'); pause; 
mypsol('T3-2','pt15');
%% soln plots, POs 
myhopl('T1-H1','pt2',[1 2]); pause; myhopl('T1-H1','pt8'); pause 
myhopl('T3-H1','pt8'); 