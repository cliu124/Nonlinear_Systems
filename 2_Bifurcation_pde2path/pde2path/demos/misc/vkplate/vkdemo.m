% script for von Karman plate
close all; clear all; 
fprintf('Init p for trivial branch of von Karman plate and\'); 
mypause('run findbif for first 3 bifurcation points'); 
p=[]; p=vkinit(p); p.sol.ds=1; p=findbif(p,3);
mypause('2 primary and one secondary bifurcations'); 
q=swibra('p','bpt1','q'); q=cont(q);
r=swibra('p','bpt2','r'); r=cont(r);
w=swibra('r','bpt1','w'); w.smod=3;w=cont(w);
mypause('Plot bifurcation diagram'); 
figure(3);clf; cmp=2;
plotbra(q,3,cmp,'ms',10,'lwst',2,'lwun',2,'cl','k');
plotbra(r,3,cmp,'ms',10,'lwst',2,'lwun',2,'cl','b');
plotbra(w,3,cmp,'ms',10,'lab',[3,6],'lw',5,'cl','r');
axis([5 8 0 8]);xlabel('\lambda'); ylabel('||u_1||_2');
mypause('Plot some solutions'); 
plotsolf('q','bpt1',4,1,2);plotsolf('r','bpt1',5,1,2);
plotsolf('w','pt3',6,1,2);plotsolf('w','pt6',7,1,2);
plotsolf('q','bpt1',8,3,2);plotsolf('r','bpt1',9,3,2);