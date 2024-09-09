% Script file for solutions plots of the tutorial "Fold and branch point
% continuation in a Schnakenberg system and details of branch plotting - a
% pde2path tutorial. The cells are named after the sections in which the
% solution plots are shown. Each cell runs independet of the others but
% requires the corresponding cmds to be run before.
% plotsol('dir','point') plots the first component of the specified point
% in figure 1, plotsol('dir','point',wnr,cmp) plots component cmp in figure
% wnr.

%% Basic fold continuation
figure(1);
clf;
plotsol('per_f1D','fpt1');
figure(2);
clf;
plotsol('per2a_f1D','fpt1',2,1);
%% Exercise 1D
figure(1);
clf;
plotsol('per1_ex','pt120');
figure(2);
clf;
plotsol('per2_ex','pt180',2,1);
figure(3);
clf;
plotsol('snak1_ex','pt300',3,1);
figure(4);
clf;
plotsol('snak2_ex','pt160',4,1);
%% 2D fold continuation
figure(1);
clf;
plotsol('hex_f2D','fpt1');
figure(2);
clf;
plotsol('hex_f2D','fpt1',2,2);