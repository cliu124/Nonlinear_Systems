%% stat. BD for sloc model: (1D, hence ly below just dummy) 
close all; keep pphome; p=[]; 
%% FSC/FSI branch; init, then cont to find bif.points 
lx=2*pi/0.44; ly=1; nx=20; sw=1; ndim=1; p=slinit(p,lx,ly,nx,sw,ndim); 
p=setfn(p,'f1'); p.nc.nsteps=10; screenlayout(p); p=cont(p,100);  
%% FSM branch 
sw=2; p=slinit(p,lx,ly,nx,sw,ndim); p=setfn(p,'f2'); 
p.nc.dsmax=0.2; p.sol.ds=0.1; p=cont(p,25); 
%% bif from f1  
p=swibra('f1','bpt1','p1',-0.05); p.nc.dsmax=0.3; p=cont(p,150); 
p=swibra('f1','bpt2','p2',-0.05); p.nc.dsmax=0.3; p=cont(p,150); 
p=swibra('f1','bpt3','p3',-0.05); p.nc.dsmax=0.3; p=cont(p,150); 
%% plotting of BD, L2-norm
figure(3); clf; pcmp=3; ms=5; 
plotbra('f1','bpt1',3,pcmp,'ms',ms,'lab',13,'cl','k','fancy',0); % FSC 
plotbra('f1','pt77',3,pcmp,'ms',ms,'fp',15,'cl','k','lab',38,'fancy',0); % FSI 
plotbra('f2','pt26',3,pcmp,'ms',ms,'lab',12,'cl','b','fancy',0); % FSM 
plotbra('p1','pt100',3,pcmp,'ms',ms,'cl','r','lab',[15 68],'fancy',0); 
plotbra('p2','pt100',3,pcmp,'ms',ms,'cl','m','lab',15,'fancy',0); 
plotbra('p3','pt100',3,pcmp,'ms',ms,'cl','g','lab',19,'fancy',0); 
 axis([0.56 0.76 0.4 1.65]);
xlabel('b'); ylabel('||v||_{2,a}');
text(0.61,1.6,'FSM','color','b','fontsize',16);text(0.56,1.59,'p1','color','r','fontsize',16);
text(0.56,1.45,'p2','color','m','fontsize',16);text(0.56,1.28,'p3','color','g','fontsize',16);
text(0.56,1,'FSI','color','k','fontsize',16);text(0.56,0.5,'FSC','color','k','fontsize',16);
%% JC
figure(3); clf; pcmp=4; ms=0; 
plotbra('f1','pt77',3,pcmp,'ms',ms,'fp',1,'cl','k','fancy',0,'lab',[13,38]); % FSI
plotbra('f2','pt26',3,pcmp,'ms',ms,'fp',1,'cl','b','fancy',0); % FSM 
plotbra('p1','pt100',3,pcmp,'ms',ms,'cl','r','lab',16,'fancy',0);
xlabel('b'); ylabel('J_{ca}'); axis([0.6 0.75 -2.7 -1.8]); 
%% solution plotting (adapt as necessary) 
plot1Df('p3','pt19',1,1,1,2); 
%% CSS-values: 
cssvalf('p1','pt74'); 

