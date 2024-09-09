%% branch plot, c0=0, E over lam1
f=3; mclf(f); xlab='\lambda_1'; c=[2 15]; ylab='E';  
plotbra('0','pt10',f,c,'cl','k'); 
plotbra('0/o',f,c,'cl',p2pc('r3'),'lab',[20]); % oblate 
plotbra('0/p',f,c,'cl',p2pc('o1'),'lsw',0); % prolate 
plotbra('0/o-1','pt27',f,c,'cl',p2pc('g1'), 'lsw',0); % oblate-D2
plotbra('0/o-2',f,c,'cl',p2pc('v2'),'lsw',0,'fms',0); % stoma!
xlabel(xlab); ylabel(ylab); box on; grid on; 
axis([3 9 0.24 0.6]); 
%% branch plot, c0=0, E over v 
f=5; mclf(f); xlab='v'; c=[18 15]; ylab='E';  
%plotbra('0','bpt1',f,c,'cl','k'); 
plotbra('0/o',f,c,'cl',p2pc('r3'),'lab',[20]); % oblate 
plotbra('0/p',f,c,'cl',p2pc('o1'),'lab',25); % prolate 
plotbra('0/o-1',f,c,'cl',p2pc('g1'), 'lab',[12, 32]); % oblate-D2
plotbra('0/o-2',f,c,'cl',p2pc('v2'),'lab',[13, 34],'fms',0); % stoma!
xlabel(xlab); ylabel(ylab); 
title('c_0=0'); box on; grid on; %'FontSize',20); %axis tight; 
axis([0.47 1 0.249 0.59]); axis([0.5 1 0.249 0.59]); 
%% branch plot, c0=0, zoom  
f=5; mclf(f); xlab='v'; c=[18 15]; ylab='E';  
plotbra('0/o',f,c,'cl',p2pc('r3'),'lab',[20]); % oblate 
plotbra('0/p',f,c,'cl',p2pc('o1'),'lsw',0); % prolate 
plotbra('0/o-1',f,c,'cl',p2pc('g1'),'lsw',0);  % oblate-D2
plotbra('0/o-2',f,c,'cl',p2pc('v2'),'lab',[13, 34],'fms',0); % stoma!
xlabel(xlab); ylabel(ylab); 
title('c_0=0'); box on; grid on; %'FontSize',20); %axis tight; 
axis([0.505 0.67 0.485 0.55]); 
%% soln plots
global p2pglob; p2pglob.vi=[60,20]; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.cm='parula'; 
p2pglob.showbd=0; p2pglob.faceal=1; p2pglob.cb=0; p2pglob.tsw=0; p2pglob.axlab=0;
plotsol('0/p','pt25'); xticks([-0.5 0.5]); yticks('auto'); pause; 
plotsol('0/o','pt20'); xticks([-0.5 0.5]); xticks('auto'); pause; 
%%
p2pglob.vi=[80,10]; 
plotsol('0/o-1','pt6'); yticks([-0.5 0.5]); pause; plotsol('0/o-1','pt32'); yticks([-0.5 0.5]);
%%
p2pglob.vi=[60,20]; p2pglob.faceal=0.25; plotsol('0/o-2','pt13');  yticks('auto'); pause; 
plotsol('0/o-2','pt34'); 
%% branch plot, c0=0, Gb over v 
f=4; mclf(f); xlab='v'; c=[18 26]; ylab='G_b';  
%plotbra('0','bpt1',f,c,'cl','k'); 
plotbra('0/o','pt32',f,c,'cl',p2pc('r3'),'lab',29); %,'lab',[29]); % oblate 
plotbra('0/p','pt20',f,c,'cl',p2pc('o1'),'lab',20); %,'lab',20); % prolate 
plotbra('0/o-1','pt35',f,c,'cl',p2pc('g1'), 'lab',[42]); % oblate-D2
plotbra('0/o-2b','pt15',f,c,'cl',p2pc('v2')); %,'lab',[25, 47],'fms',0); % stoma!
grid on; box on; xlabel(xlab); ylabel(ylab); 
%% branch plot, c0=0, dela over v 
f=4; mclf(f); xlab='v'; c=[18 21]; ylab='\Delta a';  
%plotbra('0','bpt1',f,c,'cl','k'); 
plotbra('0/o','pt32',f,c,'cl',p2pc('r3'),'lab',29); %,'lab',[29]); % oblate 
plotbra('0/p','pt20',f,c,'cl',p2pc('o1'),'lab',20); %,'lab',20); % prolate 
plotbra('0/o-1','pt35',f,c,'cl',p2pc('g1'), 'lab',[42]); % oblate-D2
plotbra('0/o-2','pt15',f,c,'cl',p2pc('v2')); %,'lab',[25, 47],'fms',0); % stoma!
grid on; box on; xlabel(xlab); ylabel(ylab); 