%% ac2D cmds2; mesh-ada via oopde and trullekrul, initialize with very coarse mesh 
par=[1 -0.2 1 0]; % parameters [c lambda gamma d]
lx=2*pi; ly=pi; sw.ref=0; sw.sym=1; p=acinit(lx,ly,22,par,sw); 
p.nc.ilam=2; p.nc.lammax=5; p.sol.ds=0.1; p.nc.dsmax=0.2; 
p=setfn(p,'trc'); p.plot.pstyle=1; p.usrlam=[2 4]; plotsol(p); 
%% find first 3 BPs via findbif
p=findbif(p,3);
%% switch to first 3 bifurcating branches and continue
p=swibra('trc','bpt1','c1',0.1); p=cont(p,30); 
p=swibra('trc','bpt2','c2',0.1); p=cont(p); 
p=swibra('trc','bpt3','c3',0.1); p=cont(p); 
%% mesh-ada-OOPDE, yielding np=1220, nt=2366
v1=[-20,30]; v2=[0,90]; 
p=loadp('c1','pt27'); p.np, plotsol(p); view(v1); pause; view(v2); p=setfn(p,'d1'); 
p=oomeshada(p,'ngen',2); view(v1);
p.pdeo.grid, pause; view(v2); fprintf('error-est: %g\n',p.sol.err); 
title(''); xlabel(''); ylabel(''); 
%% mesh-ada-trullekrul, yielding np=900, nt=1700, 
% error-est (oopde) typically very bad (if qualP=0) (recommended in 2D)
p=loadp('c1','pt27','c1-ar1'); plotsol(p); p.np
p.sw.trul=1; op=troptions2D(); op.etafac=2e-5; op.verbose=2; op.innerit=5; 
p.trop=op; p.trop.sw=15; p=oomeshada(p,'ngen',3); view(v1); pause; view(v2); 
p.file.count=27; p.fuha.savefu(p); 
[p,idx]=e2rs(p,p.u); fprintf('error-est: %g\n',p.sol.err); 
%% trullekrul, no coarsening 
p=loadp('c1','pt27','c1-ar2'); 
p.sw.trul=1; op=troptions2D(); op.etafac=2e-5; op.verbose=2; op.innerit=5; 
p.trop=op; p.trop.sw=3; p=oomeshada(p,'ngen',3); view(v1); pause; view(v2); 
p.file.count=27; p.fuha.savefu(p); 
[p,idx]=e2rs(p,p.u); fprintf('error-est: %g\n',p.sol.err); 
%% reload c1-ar2, then coarsen
p=loadp('c1-ar2','pt27','c1-ar2c'); p.trop.sw=4; 
p=oomeshada(p,'ngen',3); view(v1); pause;  view(v2); 
[p,idx]=e2rs(p,p.u); fprintf('error-est: %g\n',p.sol.err); 
%% reload c1-ar2, then coarsen and move 
p=loadp('c1-ar2','pt27','c1-ar2c'); p.trop.sw=5; 
p=oomeshada(p,'ngen',5); view(v1); pause;  view(v2); 
[p,idx]=e2rs(p,p.u); fprintf('error-est: %g\n',p.sol.err); 
%% mesh-adaption during continuation, e2rs
p=swiparf('c1','pt27','d1',4); p.sol.ds=0.1; p.nc.lammax=5; clf(2);
p.nc.amod=5; p.nc.sig=0.6; p.nc.ngen=3; p=cont(p,30); % set options for adaption and run 
%% mesh-adaption during continuation, trullekrul 
p=swiparf('c1','pt27','d1a',4); p.sol.ds=0.1; p.nc.lammax=5; p.nc.ngen=2; p.nc.amod=5; 
p.sw.trul=1; op=troptions2D(); op.verbose=2; op.etafac=2e-5; % init trullekrul-options
p.trop=op; p=cont(p,30);  % put trullekrul-options into p and run 
%% mesh-adaption during continuation, trullekrul, with additional coarsening 
p=swiparf('c1','pt27','d1b',4); p.sol.ds=0.1; p.nc.lammax=5; p.sw.trul=1; 
op=troptions2D(); op.verbose=2; op.etafac=2e-5; p.trop=op;  
op.npb=500; % desired number of points after coarsening step 
op.sw=5; p.trcop=op; % put trulle COARSENING options into p 
p.nc.ngen=2; p.nc.amod=5; p=cont(p,30); 
%% plot BD 
figure(3); clf; c=0; 
plotbra('d1','pt30',3,c,'lp',28,'fp',18,'lab',[20 23 25]); 
plotbra('d1a','pt30',3,c,'cl','r','fp',18,'lp',28,'lab',[20 23 25]);  
%% solution plots, e2rs adapatation 
ii=[23 25]; v1=[-30 30]; 
for i=ii; p=loadp('d1',['pt' mat2str(i)]); plotsol(p); 
    title(['d1/pt' mat2str(i) ', n_p=' mat2str(p.np)]); 
    view(v1); pause; 
end
%%  solution plots, trulle adapatation 
ii=[23 25]; 
for i=ii; p=loadp('d1a',['pt' mat2str(i)]); plotsol(p); 
    title(['d1a/pt' mat2str(i) ', n_p=' mat2str(p.np)]); 
    view(v1); pause; 
end  