%% demo vkplate (von Karman plate) 
% run cell-by-cell;
%%  init and findbif 3 bif.points from trivial branch
close all; %clear all; 
p=[]; p=vkinit(p); p.sol.ds=1; p=findbif(p,3);
%% 2 primary and one secondary bifurcations 
q=swibra('p','bpt1','q'); q=cont(q);
r=swibra('p','bpt2','r'); r.nc.nsteps=10; r.sw.sfem=0; r=cont(r);
w=swibra('r','bpt1','w',0.5); w.smod=3; w=cont(w,15);
%% branch plotting 
figure(3);clf; cmp=0; 
plotbra('q',3,cmp,'lwst',2,'lwun',2,'cl','k');
plotbra('r',3,cmp,'lwst',2,'lwun',2,'cl','b');
plotbra('w',3,cmp,'lab',[3,6],'lw',5,'cl','r');
axis([5 8 0 8]);xlabel('\lambda'); ylabel('||u_1||_2');
%% solution plotting 
plotsolf('q','bpt1',1,1,2);pause; plotsolf('r','bpt1',1,1,2);pause;
plotsolf('w','pt10',1,1,2);pause; plotsolf('w','bpt1',1,1,2);pause; 
plotsolf('q','bpt1',1,3,2);pause; plotsolf('r','bpt1',1,3,2);
%% some a posteriori error checks
p=loadp('q','bpt1');errcheck(p)