%% Nonlinear BC demo (nlbc) 
% $-\Delta u=0$  with b.c. $\partial_n u+\lambda s(x,y)f(u)=0$  on unit disk 
%% preparation 
close all; format compact; % clear all; 
p=[]; nx=30; p=nlbcinit(p,nx); 
%% find bifpoints from trivial branch; 
p=findbif(p,4); 
%% find bifpoints from trivial branch for negative lam 
p=setlam(p,0.1); p=setfn(p,'pb'); p=resetc(p); p.sol.ds=-0.5; p=findbif(p,4); 
%% bif. branches, 
q=swibra('p','bpt1','q1',0.02); q=cont(q,15);
q=swibra('p','bpt1','q1b',-0.02); q=cont(q,15);
q=swibra('p','bpt2','q2',0.02); q=cont(q,15);
q=swibra('pb','bpt1','r1',0.02); q=cont(q,15);
q=swibra('pb','bpt2','r2',0.02); q=cont(q,15);
q=swibra('pb','bpt3','r3',0.02); q=cont(q,15);
%% plot BD 
figure(3); clf; cmp=0; 
plotbraf('p',3,cmp,'cl','k','lwst',5); 
plotbraf('pb','pt10',3,cmp,'cl','k'); 
plotbraf('q1','pt15',3,cmp,'cl','b','lab',10); 
plotbraf('q1b','pt15',3,cmp,'cl','b','lab',5); 
plotbraf('q2','pt15',3,cmp,'cl','r','lab',10);
plotbraf('r1','pt10',3,cmp,'cl','g','lwst',2);
plotbraf('r2',3,cmp,'cl','m');
plotbraf('r3',3,cmp,'cl','c');
axis([-3 4.1 0 1.5]); xlabel('\lambda'); ylabel('||u||_2');
%% plotting solns 
a=30;b=60;ps=3; s='fontsize'; fs=16; 
plotsolf('q1','pt10',1,1,ps);xlabel('x',s,fs);ylabel('y',s,fs);zlabel(''); view(a,b);pause 
plotsolf('q1b','pt5',1,1,ps);xlabel('x',s,fs);ylabel('y',s,fs);zlabel(''); view(a,b);pause 
plotsolf('q2','pt10',1,1,ps);xlabel('x',s,fs);ylabel('y',s,fs);zlabel(''); view(a,b);

