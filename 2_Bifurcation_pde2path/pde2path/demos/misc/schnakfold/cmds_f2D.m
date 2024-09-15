%% 1 - creating a clear working space
close all; keep pphome;
%% 2 - initialising the problem (2D)
p=[];
par=[sqrt(60)*sqrt(3-sqrt(8))+1e-3, 0, 60]; % [lambda, sigma, d]
dom=[4,1]; % domain parameter
p=schnakinit(p,dom,20,par);
p.nc.dsmax=1e-1; % increase dsmax for faster calculation
p.file.smod=1; % store each solution
p=setfn(p,'tr_f2D');
% plot improvments for 2D
p.plot.pstyle=2;
p.plot.cm=hot;
p.fuha.outfu=@schnakbra_f2D; % new branch data (in particular L8-norm with ||1||=1)
kc=sqrt(sqrt(2)-1);
p.Om=16*pi^2*(dom(1)/kc*dom(2)/(sqrt(3)*kc)); % interval length
% names for cmp of schnakbra
p.plot.auxdict={'\lambda','\sigma','d','||u||_{\infty}','min(|u|)','||u||_8'};
%% 3 - find first two bif-points from homog. branch 
p.nc.nsteps=30;
p=findbif(p,2); % find first two bif points in max p.nc.nsteps steps, if possible
%% 4 - branch-switch to cold hexagons 
p=swibra('tr_f2D','bpt2','hex_f2D',0.05); % switch to cold hexagon branch
p.sw.foldcheck=1; % detect folds
p.sw.bifcheck=0; % disable bif detection
p=cont(p,10); % cont for max of 10 steps
%% 5 - fold continuation 
p=spcontini('hex_f2D','fpt1',2,'fold_f2D'); % init fold continuation in par 2
p.sol.ds=-1e-3;  % new stepsize in new primary parameter
p.plot.bpcmp=1;  % plot lam of fold position over sigma in fig 2 now 
clf(2); %p.sw.spjac=0; 
p.nc.lammin=-10; %p.nc.lammin=-0.5; % set minimal sigma (!) to -0.5
p=cont(p,5); % cont for a max of 15 steps
%% 6 - cont. in lam again from foldpoint
p=spcontexit('fold_f2D','pt10','hex2a_f2D'); % back to normal cont
p.sol.ds=-0.005;
p.nc.dsmax=0.02;
p.nc.lammin=3.2; % minimal lambda (!) is 3.2 now
p.plot.bpcmp=0; % plot l2-norm over lam in fig 2 again
clf(2);
p=cont(p,25);
p=spcontexit('fold_f2D','pt10','hex2b_f2D'); % back to normal cont
p.sol.ds=0.01;
p.nc.dsmax=0.05;
p.nc.lammin=3.208; % set min lambda to approx bif point
p.plot.bpcmp=0; % plot l2-norm over lam in fig 2 again
p=cont(p,10); % cont for a max of 10 steps  
%% 7 - plot BD 
figure(3);
clf;
% plot analytical trivial branch
plot([3.208,3.24],[3.208,3.24],'color','k','Linewidth',4);
hold on;
plot([3.208,3.05],[3.208,3.05],'color','k','Linewidth',2);
% plot computed branches with l8-norm
plotbra('tr_f2D','pt6',3,6,'cl','k');
plotbra('hex_f2D','pt9',3,6,'fp',2,'cl','b');
plotbra('hex2a_f2D','pt24',3,6,'cl','r');
plotbra('hex2b_f2D','pt9',3,6,'cl','r');
% post processing
axis([3.2 3.26 3.2 3.5]); 
text(3.202,3.33,'ch, \sigma=0','color','b','fontsize',16);
text(3.225,3.4,'ch, \sigma=-0.1395','color','r','fontsize',16);%3.205,3.432
text(3.225,3.21,'hom','color','k','fontsize',16);
%% 8 - plot lambda over sigma for fold 
figure(4);
clf;
plotbra('fold_f2D','pt16',4,1); % plot lam over sigma for fold position