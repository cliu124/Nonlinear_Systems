%% ac2Dsq; 2nd BP is double, hence use cswibra 
close all; keep pphome; 
%% init  
p=[]; par=[1 -0.2 1 0]; lx=pi; ly=lx; p=acinit(p,lx,ly,30,par); 
p.nc.neig=20; p.sw.bifcheck=2; % double Evals, use bifcheck=2
p.nc.ilam=2; p.nc.lammax=2; p.sol.ds=0.1; p.nc.dsmax=0.11; p=setfn(p,'tr');
%% cont trivial branch 
p=cont(p,30); 
%% switch to first and 3rd bifurcating branch (simple) and continue
p=swibra('tr','bpt1','b1',0.1); p.nc.dsmax=0.21; p=cont(p,20); 
p=swibra('tr','bpt3','b3',0.1); p.nc.dsmax=0.21; p=cont(p,20); 
%% cswibra at 2nd BP (double), then reset some parameters, slightly relaxed soltol
aux.soltol=1e-9; aux.ral=1; % here random initial alpha for Newton useful
p0=cswibra('tr','bpt2',aux); p0.sol.ds=0.1; p0.nc.dsmax=0.21; p0.sw.para=2; 
%% select directions and continue branches  
p=seltau(p0,2,'hs2'); p.sol.ds=-0.1; p=cont(p,20); % hor.stripes
p=seltau(p0,1,'vs2'); p.sol.ds=-0.1; p=cont(p,20); % vert.stripes
p=seltau(p0,3,'spo2'); p=cont(p,20); % spots
%% 3rd BP, simple again 
p=swibra('tr','bpt3','b3',0.1); p.nc.dsmax=0.21; p=cont(p); 
%% 4th BP, double 
p0=cswibra('tr','bpt4',aux); p0.sol.ds=0.1; p0.nc.dsmax=0.21; p0.sw.para=2; 
%% select directions and cont  
p=seltau(p0,3,'spo3');  p=cont(p,20); % spots
p=seltau(p0,2,'vs3'); p=cont(p,20); % stripes
%% 5th BP, double 
aux.soltol=1e-9; 
p0=cswibra('tr','bpt5',aux); p0.sol.ds=0.1; p0.nc.dsmax=0.21; p0.sw.para=2; 
%% select directions and cont  
p=seltau(p0,1,'5a');  p=cont(p,20); 
p=seltau(p0,5,'5c'); p=cont(p,20);
%% BD plotting, 
f=3; c=0; figure(f); clf;
plotbra('tr',f,c,'cl','k','lsw',0); plotbra('b1',f,c,'cl','b','lsw',0);
plotbra('hs2',f,c,'cl','r','lab',12); plotbra('spo2',f,c,'cl','m','lab',12); 
plotbra('b3','pt16',f,c,'cl','b','lab',9); 
plotbra('spo3',f,c,'cl','m','lab',8); 
plotbra('vs3',f,c,'cl','r','lab',9); 
plotbra('5a',f,c,'cl','m','lab',15); 
plotbra('5c',f,c,'cl','m','lab',13); 
xlabel('\lambda'); ylabel('||u||_2'); 
%% solution plots 
plotsol('spo2','pt12'); nolti; pause; plotsol('hs2','pt12'); nolti; pause; 
plotsol('vs2','pt12'); nolti; pause; plotsol('b3','pt9'); nolti; pause; 
plotsol('vs3','pt9'); nolti; pause; plotsol('spo3','pt8'); nolti; pause 
plotsol('5a','pt15'); nolti; pause; plotsol('5c','pt13'); nolti; 
%%
p=swibra('spo3','bpt1','du1'); 