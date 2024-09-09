%% branch plot, c0=2, E over v 
f=4; mclf(f); c=[18 15]; xlab='v'; ylab='E';
plotbra('2/o','pt25',f,c,'cl',p2pc('r3')); % oblate 
plotbra('2/o-1',f,c,'cl',p2pc('v1'),'lab',[25],'fms',0); % oblate-D3
plotbra('2/o-1-1a',f,c,'cl',p2pc('v3'),'lab',23); % bif from oblate-D3
plotbra('2/o-2',f,c,'cl',p2pc('r2'),'lab',[14]); % oblate-D4
%plotbra('2/o-4',f,c,'cl',p2pc('r3'),'lab',[15]); % oblate-D4
plotbra('2/p',f,c,'cl',p2pc('o1'),'lab', 16); % prolate 
plotbra('2/p-1',f,c,'cl',p2pc('g1'),'lab',[20,38]); % pear
xlabel(xlab); ylabel(ylab);  title('c_0=1.4'); 
axis([0.57 1 0.02 0.16]); grid on; box on
%% soln plots
global p2pglob; p2pglob.vi=[-5,10]; p2pglob.vi=[10,10]; 
p2pglob.edc='k'; p2pglob.cm='parula'; p2pglob.cb=0; p2pglob.faceal=1; yti=[-0.5 0.5]; 
plotsol('2/p','pt16'); yticks(yti); pause; plotsol('2/p-1','pt20'); pause, 
plotsol('2/p-1','pt38'); pause; p2pglob.vi=[10,40]; 
plotsol('2/o-1','pt25'); pause; plotsol('2/o-1-1a','pt23'); pause 
plotsol('2/o-2','pt14'); 
%%
plotsol('2/o','pt10');
%%
plotsol('2/o-3','pt6'); 
%%
p=swiparf('2/o-3','pt6','2/o3c0',[4, 2, 3 ,7,8,9 ,10,11,12]); 
p.sol.ds=0.01; p.nc.tol=1e-4; p=cont(p,2); p.nc.tol=1e-6; p=cont(p,10); 
%% only fake stable! 
p=loadp('2/o-5','pt7'); p.nc.eigref=-200; plotsol(p); 
plotspec(p,6); grid on; 
%% branch plot, c0=2, E over v 
f=4; mclf(f); c=[18 15]; xlab='v'; ylab='E'; %c=[14 15]; xlab='V'; ylab='E';
plotbra('2/o','pt55',f,c,'cl',p2pc('r3'),'bplab',1:5); % oblate 
plotbra('2/o-5','pt12',f,c,'cl',p2pc('g1'),'labi',2); % oblate 
%%
p=loadp('2/o','pt26'); p.nc.eigref=-200; p.nc.dsmax=0.1; p.nc.bisecmax=8; p=cont(p,20); 
%% cont to large c0  
p=swiparf('2/o-5','pt7','c0/o-5',[4,2,3,7,8,9,10,11,12]); p.nc.sigc=.05;
p.sol.ds=0.01; p.nc.dsmax=0.2; p.nc.tol=1e-4; p=cont(p,10); 
%%
p2pglob.vi=[-5,30]; 
plotsol('2/o','pt21'); % pause; plotsol('c0/o-5','pt10'); 