%% 1 - clear workspace
close all; keep pphome;
%% 2 - setting up problem
p=[];
par=[sqrt(60)*sqrt(3-sqrt(8))+5e-2, -0.6, 60];
p=schnakinit(p,2,100,par);
p.file.smod=10; % stores each 10th solution only
p=setfn(p,'tr_b1D');
%% 3 - contiunation of the trivial branch
p=cont(p,65); % continuation for a maximum of 65 steps
%% 4 - switch to periodic branch
p=swibra('tr_b1D','bpt1','per1_b1D',0.1); % switch to new branch with ds=0.01
p.nc.dsmax=0.1; % raise maximal stepsize
p.sw.bifcheck=0; % disables bifurcation point detection
p=cont(p,25); % continuation for a maximum of 25 steps
%% 5 - branch point cont
huclean(p); % reset figures
%p=spcontini('tr_b1D','bpt1',3,'bpcon_b1D'); % start branch point continuation in parameter 3 (d)
p=bpcontini('tr_b1D','bpt1',3,'bpcon_b1D'); p.sw.spjac=0; 
p.sw.bifcheck=0; % disable check for bifurcation
p.plot.bpcmp=1; % plot d over lambda now
p.sw.secpred=1; % use secant predictor instead of tangent! 
% continuation direction/sepsize (incl. min/max) for d
p.sol.ds=-0.1;
p.nc.dsmin=0.001;
p.nc.dsmax=1;
p.nc.lammin=0;
p.nc.lammax=100;
p=cont(p,20); % continuation for a maximum of 20 steps
%% exit BP cont 
%p=spcontexit('bpcon_b1D','pt20','per2_b1D'); % switch to continuation in lambda again
p=bpcontexit('bpcon_b1D','pt20','per2_b1D'); % switch to continuation in lambda again
%% Branch-sw. at continued BP
p=swibra('per2_b1D','bpt1','per2_b1D',-0.1); % switch to non-trivial branch
p.nc.dsmax=0.05; % reducing maximal stepsize
clf(2); p.plot.bpcmp=0; % reset ploting of bifurcation diagramm
p=cont(p,60); % continuation for a maximum of 40 steps
%% plot org. and new per. branch 
fnr=3;
figure(fnr);
clf;
plotbra('tr_b1D','pt65',fnr,0); % trivial branch
plotbra('per1_b1D','pt25',fnr,0,'cl','r'); % periodic branch for d=60
plotbra('per2_b1D','pt60',fnr,0,'cl','g'); % periodic branch for d=47.2768...
text(3.2,23.35,'d=60','color','r','fontsize',16);
text(2.6,20.9,'d=47.2768','color','g','fontsize',16);
%% plot lambda over d for branch point position
figure(4);
clf;
plotbra('bpcon_b1D',4,1);