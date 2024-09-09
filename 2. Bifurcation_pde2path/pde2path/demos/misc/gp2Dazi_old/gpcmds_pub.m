%% demo gpsol
% command templates for Gross-Pitaevski and vector GP 
% run cell-by-cell;
%% scalar GP
close all; clear all;  
p=[];p=gpinit(p); p.nc.nsteps=30; p=cont(p);
%% branch plotting 
plotbra(p,3,3,'lab',[10,15]);
%% solution plotting 
plotsolf('p','pt1',4,1,5); axis equal; % VF plot
plotsolf('p','pt10',5,1,5); axis equal; % VF plot
plotsolf('p','pt10',6,10,3); % amp plot
plotsolf('p','pt15',7,10,3); % amp plot
plotsolf('p','pt20',8,13,2); % phase plot
%% vector GP 
clear all; close all; q=[];q=vgpinit(q); q=cont(q);
%% branch plotting 
figure(3);clf;
plotbra(q,3,3,'lab',[1,10,15]);
plotbra(q,3,3,'lab',[1,10,15],'cl','r');
xlabel('\omega'); ylabel('p=max|u|/max|v| and p_1');
%% solution plotting
plotsolf('q','pt1',1,1,5); pause; % VF plot
plotsolf('q','pt10',1,1,5); % VF plot
