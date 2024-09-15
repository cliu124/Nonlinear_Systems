%% branch plot, c0=-2, E over v 
f=5; mclf(f); xlab='v'; c=[18 15]; ylab='E';  
plotbra('m2','bpt1',f,c,'cl','k'); 
plotbra('m2/o',f,c,'cl',p2pc('r3'),'lab',50); % oblate 
plotbra('m2/o-1',f,c,'cl',p2pc('v1'),'lab',30); % stoma 
plotbra('m2/p',f,c,'cl',p2pc('o1')); % prol
plotbra('m2/p-1',f,c,'cl',p2pc('g1'),'lab',[18,32]);
plotbra('m2/d',f,c,'cl',p2pc('b2'),'lab',[70],'bplab',2); %,'lp',26);
xlabel(xlab); ylabel(ylab); title('c_0=-1'); grid on; box on; 
axis([0.45 1 1 1.65]);  %axis([.45 .68 1.2 1.65]); 
%% soln plots
global p2pglob; p2pglob.vi=[30,60]; p2pglob.edc='k'; p2pglob.cm='parula'; 
p2pglob.showbd=2; p2pglob.cb=0; 
p2pglob.showN=0; p2pglob.cut=0; p2pglob.faceal=0.4; p2pglob.tsw=0; 
plotsol('m2/o','pt50'); pause; plotsol('m2/o-1','pt30'); pause 
p2pglob.vi=[-40,50]; plotsol('m2/d','pt70'); 
%%
p2pglob.faceal=1; plotsol('m2/d','bpt2'); pause 
p2pglob.vi=[-10,50]; plotsol('m2/p-1','pt18'); pause; 
plotsol('m2/p-1','pt32'); 
