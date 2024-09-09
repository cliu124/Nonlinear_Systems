%% soln plots for cmds1
plotsol('per/ini','pt1'); pause; plotsol('per/ini','pt2'); pause 
plotsol('per/stand','pt20'); pause; plotsol('per/rot0','pt20'); pause; 
plotsol('per/rot0','pt40'); pause; plotsol('per/wtstand','pt10'); 
%% branch-plotting for cmds1 (and cmds for stand2*) 
figure(3); clf; plotbra('per/stand','pt30',3,0,'cl','k','lab',20); 
plotbra('per/ini','pt2',3,0,'cl','m','lab',[1 2]); 
plotbra('per/stand','pt20',3,0,'cl',0.5*[1 1 1],'lsw',0); 
plotbra('per/rot0','pt40',3,0,'cl','b','lab',[20,40]); 
plotbra('per/wtstand','pt20',3,0,'cl','r','lab',10); 
axis([0.38,1.3,0 2.2]); %saveas(figure(3), 'figs/standbras.eps','eps');
%% plotting for cmds2
figure(3); clf; 
plotbra('per/rot','lab',25); plotsol('per/rot','pt25'); pause 
clf(3); plotbra('per/wtmove','lsw',0); pause 
clf(3); plotbra('per/wtasym','lab',20); plotsol('per/wtasym','pt20');
%% plotting for cmds3 - NBC 
figure(3); clf; plotbra('nbc/stand','pt20',3,0,'cl','black','lab',20); pause
plotsol('nbc/stand','pt20'); pause; plotsol('nbc/rot','pt10'); pause; 
plotsol('nbc/dom','pt10'); pause; 
plotsol('nbc/move-','pt20'); pause; plotsol('nbc/move+','pt20'); 
%% BD plotting 
clf(3);plot([-0.03 -0.02],[0 0],'-k',[-0.03 -0.02],[0 0],'-b');  legend('s','\nu'); 
plotbra('nbc/move-','pt20',3,2,'cl','k','lab',20,'lw',2);
plotbra('nbc/move-','pt20',3,3,'cl','b','lab',20,'lw',2); 
plotbra('nbc/move+','pt20',3,2,'cl','k','lab',20,'lw',2); 
plotbra('nbc/move+','pt20',3,3,'cl','b','lab',20,'lw',2); 
axis([-0.02 0.02 -0.011 0.011]); ylabel(''); 