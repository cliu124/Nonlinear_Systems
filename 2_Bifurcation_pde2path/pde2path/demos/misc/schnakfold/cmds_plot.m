%% 1 - creating a clear working space
close all; keep pphome;
%% 2 - initialising the problem
p=[];
par=[sqrt(60)*sqrt(3-sqrt(8))+5e-2, -0.6, 60]; % [lambda, sigma, d]
p=schnakinit(p,4,300,par);
p.plot.pmod=10; % shows each 10th solution in fig 2 only
p.file.smod=10; % stores each 10th solution only
p=setfn(p,'tr_plot');
%% 3 - contiunation of the trivial branch
p=cont(p,20); % continuation for a maximum of 20 steps
%% basic/cell 1 - switch to periodic branch
p=swibra('tr_plot','bpt1','per_plot',1e-2); % switch to new branch with ds=0.01
p=cont(p,150); % continuation for a maximum of 150 steps
%% basic 2 - basic plot
figure(3);
clf;
plotbra('tr_plot'); % trivial branch
% plotbra('tr_plot','pt20'); % plots the same, 'pt10' plots a shortend branch
plotbra('per_plot'); % periodic branch
% plotbra(p); % generates the same as plotbra('per_plot')
%% cell 2 - plot with options through a cell array
figure(3);
clf;
basic={'ms',2,'fancy',0,'lsw',15};
plotbra('tr_plot',basic); % trivial branch
% plotbra('tr_plot','ms',2,'fancy',0,'lsw',15); % plots the same
plotbra('per_plot',basic); % periodic branch
% plotbra(p,basic) % plots the same as plotbra('per_plot',basic)
%% usrlam 1 - usrlam data generation
p=swibra('tr_plot','bpt1','per_plot_usrlam',1e-2); % switch to periodic branch
p.usrlam=3:0.1:3.5; % set userlambdas
p=cont(p,200); % continuation for a maximum of 200 steps
%% usrlam 2 - usrlam plotting
figure(3);
clf;
plotbra(p); % periodic branch with usrlabels, as lsw=1 by default
% plotbra('per_plot'); would plot the same
plotbra('tr_plot'); % trivial branch - no labels, as p.usrlam={} in Cell 3