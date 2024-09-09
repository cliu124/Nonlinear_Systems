%% BD for bd-harvesting of LV fish; find branches by initial guesses
close all; keep pphome;  
%% rem: par=[beta b2 b3 d1 d2 ga1 ga2 rho al1 al2 p1 p2 c1 c2]; 
%            1     2  3  4  5  6   7   8   9  10  11 12 13 14
p=[]; lx=20; nx=50; sw=1; p=lvinit(p,lx,nx,sw); p0=p;
%% cont in c1
p=p0; p=setfn(p,'c1'); p.usrlam=[0.5 1 5 10]; 
p.nc.tol=1e-8; p.nc.ilam=13; p.sol.ds=0.05; p.sol.restart=1; p.nc.lammax=20;
p.sw.bifcheck=0; p.file.smod=5; p.nc.dsmax=4; p=cont(p,30);
%% plotting of BD, [J_c,J1,J2, v1,v2, l1,l2,k1,k2, h1,h2,p1-ga1*l1; p2-ga2*l2, min,max]
figure(3); clf; ms=5; 
plot([0 0.1], [-1 -1],'b'); hold on; plot([0 0.1], [-1 -1],'k'); plot([0 0.1], [-1 -1],'r'); 
leg=legend('J','J_1', 'J_2'); 
plotbra('c1','pt31',3,1,'ms',ms,'lab',[0,31],'cl','b');  
plotbra('c1','pt31',3,2,'ms',ms,'cl','k','lsw',0);
plotbra('c1','pt31',3,3,'ms',ms,'cl','r','lsw',0); 
xlabel('c_1'); ylabel('');
%% lambdas
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('\lambda_1', '\lambda_2'); 
plotbra('c1','pt31',3,6,'ms',ms,'cl','k');  
plotbra('c1','pt31',3,7,'ms',ms,'cl','r');  
xlabel('c_1'); ylabel(''); 
%% v 
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('v_1','v_2'); plotbra('c1','pt31',3,4,'ms',ms,'cl','k');  
plotbra('c1','pt31',3,5,'ms',ms,'cl','r');  
xlabel('c_1'); axis([0.1 10 0 1.4]); ylabel('');
%% k 
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('k_1','k_2'); plotbra('c1','pt31',3,8,'ms',ms,'cl','k');  
plotbra('c1','pt31',3,9,'ms',ms,'cl','r');  
xlabel('c_1'); axis([0.1 10 0 100]); ylabel('');
%% h
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('h_1','h_2'); plotbra('c1','pt31',3,10,'ms',ms,'cl','k');  
plotbra('c1','pt31',3,11,'ms',ms,'cl','r');  
xlabel('c_1'); axis([0.1 10 0 8]); ylabel('');
%% soln plot 
plotsolf('c1','pt0',1,1,10); pause; plotsolf('c1','pt31',1,1,10); 
%% ****************************** c2 *************************************
p=p0; p=setfn(p,'c2'); p.usrlam=[0.5 1 5 10]; 
p.nc.tol=1e-8; p.nc.ilam=14; p.sol.ds=0.05; p.sol.restart=1; p.nc.lammax=20;
p.sw.bifcheck=0; p.file.smod=5; p.nc.dsmax=4; p=cont(p,30);
%% plotting of BD, [J_c,J1,J2, v1,v2, l1,l2,k1,k2, h1,h2,p1-ga1*l1; p2-ga2*l2, min,max]
figure(3); clf; pcmp=1; ms=5; 
plot([0 0.1], [-1 -1],'b'); hold on; plot([0 0.1], [-1 -1],'k'); plot([0 0.1], [-1 -1],'r'); 
leg=legend('J','J_1', 'J_2'); 
plotbra('c2','pt31',3,1,'ms',ms,'lab',[0,10,31],'cl','b');  
plotbra('c2','pt31',3,2,'ms',ms,'cl','k');
plotbra('c2','pt31',3,3,'ms',ms,'cl','r');  
xlabel('c_2'); ylabel('');
%% lambdas
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('\lambda_1', '\lambda_2'); 
plotbra('c2','pt31',3,6,'ms',ms,'cl','k');  
plotbra('c2','pt31',3,7,'ms',ms,'cl','r');  
xlabel('c_2'); ylabel(''); 
%% v 
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('v_1','v_2'); plotbra('c2','pt31',3,4,'ms',ms,'cl','k');  
plotbra('c2','pt31',3,5,'ms',ms,'cl','r');  
xlabel('c_2'); axis([0.1 10 0 0.3]); ylabel('');
%% k 
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('k_1','k_2'); plotbra('c2','pt31',3,8,'ms',ms,'cl','k');  
plotbra('c2','pt31',3,9,'ms',ms,'cl','r');  
xlabel('c_2'); axis([0.1 10 0 100]); ylabel('');
%% h
figure(3); clf; plot([0 0.1], [-1 -1],'k'); hold on; plot([0 0.1], [-1 -1],'r'); 
leg=legend('h_1','h_2'); plotbra('c2','pt31',3,10,'ms',ms,'cl','k');  
plotbra('c2','pt31',3,11,'ms',ms,'cl','r');  
xlabel('c_2'); axis([0.1 10 0 8]); ylabel('');
%% soln plot 
plotsolf('c2','pt31',1,1,10); pause; plotsolf('c2','pt31',1,1,10); 