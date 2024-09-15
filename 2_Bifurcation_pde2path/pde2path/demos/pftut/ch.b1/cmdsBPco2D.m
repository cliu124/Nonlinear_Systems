%% extra commands 1D, BP cont 
close all; keep pphome; p=[]; 
%%  continuation in eps of BPT2 (simple) on 2D: 
p=[]; figure(2); clf; p=bpcontini('2D','bpt2',2,'bpc2'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.lammin=0.01; 
p.nc.del=1e-2; p.sw.spjac=1; p.fuha.spjac=@bpjac; 
p.nc.tol=1e-6; p.file.smod=1; tic; p=cont(p,10); toc
%% continue in other direction, natural fold at (eps,m)=(0.31,0)
p=loadp('bpc2','pt0','bpc2-b'); p.nc.dsmax=0.05; p.sol.ds=-p.sol.ds; p.file.smod=10; p=cont(p,20); 
%% branch plot 
mclf(3); plotbra('bpc2','pt4',3,1,'cl','b'); xlabel('\epsilon'); 
%% switch back to regular continuation
p=bpcontexit('bpc2','pt4','pbpc2'); p=resetc(p); p.nc.dsmax=0.1; 
%p.plot.bpcmp=5; p.nc.lammin=-5; p.sol.ds=1e-3; p.sw.spcalc=1; clf(2); huclean(p); p=cont(p,2); 
%% check that BP is fine, i.e., that branch-switching works
p=swibra('pbpc2','bpt1','2D2-smalleps',0.1); pause; 
p.nc.lammin=-5; p.plot.bpcmp=5; mclf(2); p=cont(p,40); 
%%  continuation in eps of BPT1 (double) on 2D: also works, but choice of kernel vector 
% is somewhat random, i.e., depends on the discretization! 
p=[]; figure(2); clf; p=bpcontini('2D','bpt1',2,'bpc3'); p.sol.ds=-0.01; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.lammin=0.01; 
p.nc.del=1e-2; p.sw.spjac=1; p.fuha.spjac=@bpjac; 
p.sw.spqjac=1; p.fuha.spqjac=@stanspqjac; p.nc.tol=1e-7; p.file.smod=1;
tic; p=cont(p,10); toc
%% switch back to regular continuation
p=bpcontexit('bpc3','pt4','pbpc3'); p=resetc(p); p.nc.dsmax=0.1; 
%% check that BP is fine, i.e., that branch-switching works
p=swibra('pbpc3','bpt1','2D1-smalleps',0.1); p.nc.lammin=-5; p.plot.bpcmp=5; mclf(2); p=cont(p,40); 
%% compare old and new branches, smaller eps give better approximations of interface lengths 
figure(3); clf; f=3; c=5; plotbra('2D2','pt40',f,c,'cl','m','lab',50); 
plotbra('2D2-smalleps','pt40',f,c,'cl',p2pc('g1')); 
plotbra('2D1-sp','pt40',f,c,'cl','b'); 
plotbra('2D1-smalleps','pt40',f,c,'cl',p2pc('b2')); 
ylabel('E_\epsilon'); 