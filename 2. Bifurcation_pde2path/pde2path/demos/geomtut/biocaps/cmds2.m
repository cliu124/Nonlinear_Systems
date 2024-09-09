%% experiments with continuation of non-axi branches to other b and lam1
% uncomment as desired
global p2pglob; p2pglob.edc='k'; p2pglob.showbd=2; p2pglob.vi=[70,30]; 
p=swiparf('c01b-1q','pt9','t2/bc1',[4 5]); 
p.nc.tol=1e-6; p.sol.ds=0.025; p=cont(p,20); 
%%
p=swiparf('c01b-1q','pt16','t2/bc2',[4 5]); p.nc.lammin=-4.5; p.nc.bisecmax=4; 
p.nc.tol=1e-6; p.sol.ds=0.025; p=cont(p,20); 
%% cont in lam1
p=swiparf('c01b-1q','pt9','t2/lc1',[2 5]); p.nc.tol=1e-6; p.sol.ds=0.01; p=cont(p,10); 
%%
p=swiparf('c01b-1q','pt16','t2/lc2',[2 5]); p.nc.tol=1e-5; p.sol.ds=0.0025; p=cont(p,15); 
%% branch plot b-cont of non-axi 
f=4; mclf(f); xlab='b'; c=8; ylab='E'; %c=7; ylab='A'; c=10; ylab='K'
plotbra('t1/bc1','pt16',f,c,'cl','k','lab',[0 10 16],'lp',16); %[15 30 55]); 
plotbra('t1/bc2',f,c,'cl','b','lab',[0 8 14]);
xlabel(xlab); ylabel(ylab); grid on
%% soln plots
p2pglob.edc='none'; p2pglob.tsw=0; p2pglob.vi=[60 10]; 
pplot('t1/bc1','pt0'); zticks([0 0.4 0.8]); pause; pplot('t1/bc1','pt10'); zticks([0 0.2]); pause; 
pplot('t1/bc1','pt16'); zticks('auto'); pause; 
%%
pplot('t1/bc2','pt0'); pause; pplot('t1/bc2','pt8'); pause; pplot('t1/bc2','pt14'); 
%% branch plot lam1-cont of non-axi 
f=5; mclf(f); xlab='\lambda_1'; c=8; ylab='E'; %c=7; ylab='A'; c=11; ylab='\delta_{mesh}'; 
plotbra('t1/lc1',f,c,'cl','k','lab',[0 8],'lp',8); 
plotbra('t1/lc2',f,c,'cl','b','lab',[0 11 15]);
xlabel(xlab); ylabel(ylab); grid on
%% soln plots
pplot('t1/lc1','pt8'); pause; 
pplot('t1/lc2','pt11'); zticks([-0.4 0 0.4 0.8]); pause; pplot('t1/lc2','pt15'); zticks('auto');
