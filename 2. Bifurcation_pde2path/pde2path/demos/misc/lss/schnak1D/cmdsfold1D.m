%% 1 - creating a clear working space
close all; keep pphome;
%% 2 - initialising the problem
p=[];
par=[sqrt(60)*sqrt(3-sqrt(8))+5e-2, -0.6, 60]; % [lambda, sigma, d]
pp=20; np=5000; p=schnakinit(p,pp,np,par);
p.plot.pmod=10; % shows each 10th solution in fig 2 only
p.file.smod=10; % stores each 10th solution only
%% 3 - contiunation of the trivial branch
tic; p=cont(p,20); toc % continuation for a maximum of 20 steps
%% 4 - switch to periodic branch and continuation with fold detection
p=swibra('p','bpt1','periodic',1e-2); % switch to new branch with ds=0.01
%p=setbel(p,0,1e-3,5,@lss); 
p=setbelilup(p,0,1e-3,5,1e-5,50); 
p.sw.verb=2; p.sw.spcalc=0; p.plot.pmod=50; 
p.sw.foldcheck=0; % enables detection of folds
p.sw.bifcheck=0; % disable detection of bifs
tic; p=cont(p,50); toc % continuation for a maximum of 150 steps
%% 5 - fold continuation in sigma
p=spcontini('periodic','fpt1',2,'fold'); % switch to fold cont in parameter sigma
p.sol.ds=-1e-3; % continue backward in sigma
clf(2);
p.plot.bpcmp=1; % plot lambda position of fold over sigma in fig 2 now
p=cont(p,50); % continuation for a maximum of 50 steps
%% 6 - continuation at new fold point in lam again
p=spcontexit('fold','pt51','periodic2a'); % exit fold continuation
p.sol.ds=1e-2; % continue forward in lambda
clf(2);
p.plot.bpcmp=0; % plot L2-norm over lambda in fig 2
p=cont(p,20); % cont for a maximum of 20 steps
p=spcontexit('fold','pt51','periodic2b'); % exit fold continuation
p.sol.ds=-1e-2; % continue forward in lambda
p.nc.lammin=3.2084; % minimal lambda - approx. bif point from trivial branch
p.plot.bpcmp=0; % plot L2-norm over lambda in fig 2
p=cont(p,120); % cont for a maximum of 120 steps
%% 7 - plot BD
figure(3);
clf;
plotbra('p'); % trivial branch
plotbra('periodic'); % periodic nranch for sigma=-0.6
plotbra('periodic2a','cl','r'); % periodic branch for sigma=-0.7617
plotbra('periodic2b','cl','r','auxdict',{'\lambda','\sigma','d','||u_1||_{\infty}','min(|u_1|)'}); % periodic branch for sigma=-0.7617 and labels axis
%% 8 - plot lambda over sigma for fold 
figure(4);
clf;
plotbra('fold',4,1,'auxdict',{'\lambda','\sigma','d','||u_1||_{\infty}','min(|u_1|)'}); % plot lambda position of fold over sigma and labels axis