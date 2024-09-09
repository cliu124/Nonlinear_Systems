%% demo chemotax, quasilinear system, uses getKuvd to get \pa_u(div(c(u)grad v))
%% Init
close all; lx=0.5; ly=2; nx=20; ny=round(ly*nx/lx);  nx=20;
p=[]; p=cheminit(p,lx,ly,nx,ny); p.sw.bifcheck=2; p.nc.mu2=0.1; p.np, pause 
p=setfn(p,'tr'); p.sw.spcalc=1; p.sw.verb=2; 
p.nc.bisecmax=5; p.nc.dsmax=0.5;  
p.nc.nsteps=100; p.sw.jac=0; p=findbif(p,8); p=cont(p,20); 
%% BP1, 1 stripe, switching off bifdetec for speed, 
p=swibra('tr','bpt1','q1',0.25); p.sw.jac=1; 
p.sw.bifcheck=0; p.sw.spcalc=0; p.nc.dsmax=0.5; p=cont(p,20); 
%% 2 stripes 
p=swibra('tr','bpt2','q2',0.5); p.sw.jac=1; 
p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p,35);
%% BP3, multiple, only compute kernel, then use gentau 
aux=[]; aux.besw=0; aux.m=3; p0=cswibra('tr','bpt3',aux); 
p0.sw.bifcheck=0; p0.sw.spcalc=0; p0.sw.jac=1; 
%% 3 stripes 
p=gentau(p0,[0 1],'q3'); p.sol.ds=0.05; p=cont(p,20);
%% single vertical stripe 
p=gentau(p0,[0 0 1],'q4'); p.sol.ds=0.05; p=cont(p,30);
%% BP6, 2 spots with a  loop along the way 
p=swibra('tr','bpt6','q5',0.05); p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p,40); 
%% 3 spots 
p=swibra('tr','bpt7','q6',0.05); p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p,25); 
%% BP7, multiple, only compute kernel, subsequently use gentau
aux=[]; aux.besw=0; aux.m=2; p0=cswibra('tr','bpt8',aux);
p0.sw.bifcheck=0; p0.sw.spcalc=0; p0.sw.jac=1; 
%% 4 spots
p=gentau(p0,[1 0],'q7'); p.sol.ds=0.05; p=cont(p,25);
%% 3 stripes
p=gentau(p0,[0 1],'q8'); p.sol.ds=0.05; p=cont(p,25);
%% plot bifurcation diagram
figure(4);clf(4);cmp=1;
plotbra('tr',4,cmp,'cl','k','ms',0); 
plotbra('q1',4,cmp,'cl',p2pc('b1'),'lab',5,'lp',20);
plotbra('q2',4,cmp,'cl',p2pc('b2'),'lab',10); 
plotbra('q3',4,cmp,'cl','b'); %,'lab',20);
plotbra('q4',4,cmp,'cl',p2pc('b3'),'lab',15); 
plotbra('q5',4,cmp,'cl',p2pc('r2'),'lab',30); 
plotbra('q6',4,cmp,'cl',p2pc('r3'),'lab',15);
plotbra('q7',4,cmp,'cl',p2pc('o1'),'lab',10);
plotbra('q8',4,cmp,'cl',p2pc('b3')); %,'lab',10);
axis([10 25 0 0.8]);xlabel('\lambda');
ylabel('||u_1-1||_{L^1}/|\Omega|');
%% plot some solns 
psol('q1','pt5',1,1,2);pause; psol('q2','pt10',1,1,2); pause 
psol('q3','pt15',1,1,2);pause; 
psol('q4','pt15',1,1,2);pause; 
psol('q5','pt30',1,1,2); pause; psol('q6','pt15',1,1,2); pause
psol('q7','pt10',1,1,2); %pause; psol('q8','pt10',1,1,2); 
%% check the diff matrices 
n=p.np; plotsol(p,1,1,1); ux=p.u; 
ux(1:n)=p.mat.Dy*p.u(1:n); %ux(1:n)=smu(p,ux(1:n)); %ux(1:n)=p.mat.M(1:n,1:n)\ux(1:n); 
plotsolu(p,ux,3,1,1); 
%% check Jacobians
p.nc.del=1e-6; 
[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); figure(4); clf; spy(Gd>e1/2);