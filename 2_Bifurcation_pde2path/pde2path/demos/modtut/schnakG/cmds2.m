%% Schnak on larger SFN (scale free network)  
close all; keep pphome; 
%% init 
p=[]; par=[10,0,0.4,200]; 
sw=4; np=250; sw=0; np='G2'; 
p=schnakinitG(p,par,np,sw); p=setfn(p,'tr2'); p.nc.lammax=12; p.nc.dsmax=0.25; 
p.sw.bifcheck=2; p.sw.verb=2; p.nc.mu1=2; p.nc.mu2=0.5; p.sol.ds=-0.1; p.nc.ilam=1; 
p.fuha.outfu=@hobra; p.plot.bpcmp=6; p.vol=1; plotsol(p); p.dfac=0.2; 
p=cont(p,20); 
%% swibra;  
p=swibra('tr2','bpt1','b1',0.5); pause; p.sw.usrlam=[]; 
p.file.smod=50; p.sw.foldcheck=0; p.sw.bifcheck=0; p=cont(p,500);  
%G=p.G; save('G2','G'); % if you like the graph, save it here  
%% swibra
p=swibra('tr2','bpt2','b2',0.1); pause;  p.sw.usrlam=[]; 
p.file.smod=50; p.sw.foldcheck=0; p.sw.bifcheck=0; p=cont(p,200); 
%% BD plot 
f=3; c=6; figure(f); clf;
plotbra('tr2','pt20',f,c,'cl','k','lsw',0); 
plotbra('b1','pt400',f,c,'cl','b','lab',[50 200 400]); 
plotbra('b2',f,c,'cl','r','lab',[100 200]); ylabel('max u'); 
%% soln plots 
plotsol('b1','pt50'); pause; plotsol('b1','pt200'); pause; plotsol('b1','pt400');
%%
plotsol('b2','pt50'); pause; plotsol('b2','pt100'); pause; plotsol('b2','pt200');
