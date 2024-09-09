%% extra commands 1D, FPT and BPT continuation (with nq=1)
close all; keep pphome; p=[]; 
%%  continuation in eps of BPT1 on 1Dhom: 
p=[]; mclf(2); p=bpcontini('1D','bpt1',2,'bpc1'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.lammin=0.01; 
p.nc.del=1e-3; p.sw.spjac=1; p.fuha.spjac=@bpjac; p.fuha.outfu=@stanbra; p.file.smod=1; 
p.nc.tol=1e-6; tic; p=cont(p,10); toc
%% continue in other direction, natural fold at (eps,m)=(0.31,0)
p=loadp('bpc1','pt0','bpc1-b'); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p.file.smod=10; p=cont(p,80); 
%% branch plot 
mclf(3); plotbra('bpc1','pt8',3,1,'cl','b','lab',7); xlabel('\epsilon'); box on; grid on
%% switch back to regular continuation
rp='pt7'; od='pbpc1'; p=bpcontexit('bpc1',rp,[od '-a']); p=resetc(p); p.nc.dsmax=0.1; 
p.plot.bpcmp=5; p.nc.lammin=-5; p.sol.ds=1e-3; p.sw.spcalc=1; clf(2); huclean(p); 
p.fuha.outfu=@Jbra; p=cont(p,2); 
%% check that BP is fine, i.e., that branch-switching works
p=swibra('pbpc1-a','bpt1','a1-1',0.1); p.nc.lammin=-1; p.plot.bpcmp=5; p.fuha.outfu=@Jbra; 
p=cont(p,40); 
%% plot new branches; E_\eps poor due to coarse resolution! 
f=4; c=5; mclf(f); 
%plotbra('1D','pt30',f,c,'cl','k','lsw',0); 
plotbra('a1','pt50',f,c,'cl','b','lsw',0); 
plotbra('a1-1','pt40',f,c,'cl','m','lsw',0); 
ylabel('E_\epsilon'); axis([-0.95 -0.4 0.93 1.1]); 
%%
plotsol('a1-1','fpt1'); 
plotsol('a1-1','fpt1'); 