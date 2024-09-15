%% 1 - creating a clear working space
close all; keep pphome;
%% 2 - initialising the problem
p=[];
par=[sqrt(60)*sqrt(3-sqrt(8))+5e-2, -0.6, 60]; % [lambda, sigma, d]
p=schnakinit(p,4,300,par);
p.plot.pmod=10; % shows each 10th solution in fig 2 only
p.file.smod=10; % stores each 10th solution only
p=setfn(p,'tr_f1D');
%% 3 - contiunation of the trivial branch
p=cont(p,20); % continuation for a maximum of 20 steps
%% 4 - switch to periodic branch and continuation with fold detection
p=swibra('tr_f1D','bpt1','per_f1D',1e-2); % switch to new branch with ds=0.01
p.sw.foldcheck=1; % enables detection of folds
p.sw.bifcheck=0; % disable detection of bifs
p=cont(p,150); % continuation for a maximum of 150 steps
%% 5 - fold continuation in sigma
p=spcontini('per_f1D','fpt1',2,'fold_f1D'); % switch to fold cont in par. sigma
p.sol.ds=-1e-3; % continue backward in sigma
figure(2); clf; %p.sw.spjac=0; 
p.plot.bpcmp=1; % plot lambda position of fold over sigma in fig 2 now
p=cont(p,50); % continuation for a maximum of 50 steps
%% 6 - continuation at new fold point in lam again
p=spcontexit('fold_f1D','pt100','per2a_f1D'); % exit fold continuation
p.sol.ds=1e-2; % continue forward in lambda
figure(2); clf; p.plot.bpcmp=0; % plot L2-norm over lambda in fig 2
p=cont(p,20); % cont for a maximum of 20 steps
%%
p=spcontexit('fold_f1D','pt100','per2b_f1D'); % exit fold continuation
p.sol.ds=-1e-2; % continue backward in lambda
%p.nc.lammin=3.2084; % minimal lambda - approx. bif point from trivial branch
p.plot.bpcmp=0; % plot L2-norm over lambda in fig 2
p=cont(p,120); % cont for a maximum of 120 steps
%% 7 - plot BD
figure(3);
clf;
plotbra('tr_f1D'); % trivial branch
plotbra('per_f1D'); % periodic nranch for sigma=-0.6
plotbra('per2a_f1D','cl','r'); % periodic branch for sigma=-0.7617
plotbra('per2b_f1D','cl','r'); % periodic branch for sigma=-0.7617
%% 8 - plot lambda over sigma for fold 
figure(4);
clf;
plotbra('fold_f1D',4,1); % plot lambda position of fold over sigma