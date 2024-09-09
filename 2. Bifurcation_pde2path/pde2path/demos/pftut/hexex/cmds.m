%% hexex cmds 
close all; keep pphome; 
%% trivial branch, with artificial increase of lam to near BPs 
p=[]; lam=6; p=hexinit(p,0.05,lam);  p=setfn(p,'tr'); 
p=cont(p,10); p.u(p.nu+1)=17; p=cont(p,10); p.u(p.nu+1)=30; p=cont(p,30); 
%% cswibra detects non-isol. zeros 
aux=[]; p0=cswibra('tr','bpt2',aux); 
%% use gentau to generate branches
p0.pm.resfac=1e-3; p0.nc.dsmin=0.1; p0.nc.dsmax=1; p0.sol.ds=0.1; 
p=gentau(p0,[1 0],'b1'); p=pmcont(p,30); 
%p=gentau(p0,[0 1],'b2'); p=pmcont(p,30); 
%% mixed modes fall to one of the isotropy classes 
p=gentau(p0,[1 2],'b3');  p.nc.tol=1e-6; p=cont(p,10); 
%% plot BD (though L2-norms almost identical) 
figure(3); clf; pcmp=0; 
plotbra('b1','pt30',3,pcmp,'cl','r','lab',30); 
plotbra('b2','pt30',3,pcmp,'cl','b','lab',30); 
%% solution plots 
plotsol('b1','pt30',21,1,2); pause; plotsol('b2','pt30',22,1,2);
%% check of meshada
p=loadp('b1','pt30'); p.nc.sig=0.05; p=oomeshada(p); 