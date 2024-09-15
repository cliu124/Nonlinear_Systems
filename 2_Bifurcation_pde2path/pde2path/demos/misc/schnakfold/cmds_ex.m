%% 1 - creating a clear working space
close all; keep pphome;
%% 2 - initialising the problem
p=[];
par=[sqrt(60)*sqrt(3-sqrt(8))+5e-2, -0.6, 60]; % [lambda, sigma, d]
p=schnakinit(p,4,300,par);
p.nc.lammin=3;
p.nc.lammax=3.5; %min/max lambda
p.file.smod=20; % save each 20 point only
p.plot.pmod=20; % plot each 20 solution only
p=setfn(p,'tr_ex'); % set problem directory to p_ex
%% 3 - continuation of the trivial branch
p=cont(p,50); % continuation for a maximum of 50 steps
p=loadp('tr_ex','pt0','tr2_ex'); % load starting point
p.sol.ds=-p.sol.ds; % change direction of continuation
p=cont(p,50); % contiunation for max 50 steps
%% 4 - continue the periodic branch with wavenumber kc
p=swibra('tr_ex','bpt1','per1_ex'); % switch to periodic branch with wavenumber kc
p=cont(p,500); % cont for max 500 steps
%% 5 - continue the first snaking branch from branch with wavenumber kc
p=swibra('per1_ex','bpt1','snak1_ex'); % switch to localised pattern branch
p.sw.foldcheck=1; % needed for fold-cont only
p=cont(p,700); % cont for max 700 steps
%% BD 6 - continue a periodic branch with other wavenumber
p=swibra('tr_ex','bpt2','per2_ex'); % switch to periodic branch with wavenumber <kc
p=cont(p,500); % cont for max 500 steps
%% BD 7 - continue the second snaking branch of kc branch
p=swibra('per1_ex','bpt2','snak2_ex'); % switch to another localised pattern branch
%p=setbel(p,0,1e-3,10,@lss); % this may speed up things
tic; p=cont(p,400); toc % cont for max 400 steps
%% BD 8 - plot of bifurcation diagramm
figure(3);
clf;
plotbra('tr_ex','bplab',[1,2],'fancy',0); % plot trivial branch
plotbra('tr2_ex','fancy',0); % plot trivial branch
plotbra('per1_ex','bplab',[1,2],'cl','r'); % plot periodic branch wavenumber kc
plotbra('per2_ex','bplab',[6,7],'cl','m'); % plot periodic branch wavenumber <kc
plotbra('snak1_ex','bpt1','cl','b','fms',0); % plot localised pattern branch
plotbra('snak2_ex','bpt5','cl','c'); % plot second localised pattern branch
%% FC 6 - fold contiunation of the right snaking side
clf(2);
p=spcontini('snak1_ex','fpt3',2,'fr1_ex'); % start fold continuation in sigma at fpt3
p.sol.ds=-1e-2; % negativ direction in sigma cont
p.nc.lammin=-0.85;
p.plot.bpcmp=1; % plot lambda of fold over sigma now
p=cont(p,500); % cont for a max of 500 steps
p=spcontini('snak1_ex','fpt3',2,'fr2_ex'); % start fold continuation in sigma at fpt3
p.sol.ds=1e-2; % positiv direction in sigma cont
p.nc.lammin=-0.85;
p.nc.lammax=-0.5;
p.plot.bpcmp=1; % plot lambda of fold over sigma now
p=cont(p,500); % cont for a max of 500 steps
%% FC 7 - fold contiunation of the left snaking side
p=spcontini('snak1_ex','fpt4',2,'fl1_ex'); % start fold continuation in sigma at fpt4
p.sol.ds=-1e-2; % negativ direction in sigma cont
p.nc.lammin=-0.85;
p.plot.bpcmp=1; % plot lambda of fold over sigma now
p=cont(p,500); % cont for a max of 500 steps
p=spcontini('snak1_ex','fpt4',2,'fl2_ex'); % start fold continuation in sigma at fpt4
p.sol.ds=1e-2; % positiv direction in sigma cont
p.nc.lammin=-0.85;
p.nc.lammax=-0.5;
p.plot.bpcmp=1; % plot lambda of fold over sigma now
p=cont(p,500); % cont for a max of 500 steps
%% FC 8 - plot fold continuation lambda over sigma for fpt3 (black) and fpt4(red)
figure(10);
clf;
plotbra('fr1_ex',10,1); % right fold continuation
plotbra('fr2_ex',10,1); % right fold continuation
plotbra('fl1_ex',10,1,'cl','r'); % left fold continuation
plotbra('fl2_ex',10,1,'cl','r'); % left fold continuation