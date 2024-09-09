close all; keep pphome; 
%% C1: brusselator on small 2D domain, init, and cont across Hopf line 
p=[]; lx=sqrt(2)*pi/0.7; ly=1*lx; nx=11; a=3; b=8; du=6; dv=10; al=0; om=1; 
par=[a b du dv al om]; sw.sym=1; p=bruinit(p,[lx ly],nx,par,sw); 
p=initeig(p,10); p.nc.neig=[4, 4]; p.file.smod=5; p=setfn(p,'hom2Db'); 
p.nc.ilam=2; p.sol.ds=0.2; p=findbif(p,2); p=cont(p,10);  
%% C2: Hopf bif (hom), first with hopf.bisec=2 for speed 
figure(2); clf; aux.tl=30; p=hoswibra('hom2Db','hpt1',0.1,4,'2dh1',aux); 
p.hopf.bisec=2; p.nc.dsmax=0.3; p.nc.tol=1e-6; 
p.hopf.fltol=1e-2; bw=1; beltol=1e-6; belimax=5; 
p=setbel(p,bw,beltol,belimax,@lss); p=cont(p,40);
%% C3: reload and find first 2 PDs using a smaller dsmax and larger hopf.bisec   
p=loadp('2dh1','pt10'); p.sol.ds=0.05; p.nc.dsmax=0.05; p.hopf.bisec=5; p=cont(p,10); 
%% C4: PD bif to primary oscillating squares
ds=0.5; aux=[]; aux.sw=-1; p=poswibra('2dh1','bpt1','2dpd1',ds,aux); 
p.nc.dsmax=0.5; p.hopf.bisec=2; p.nc.tol=1e-6; p=cont(p,30); 
%% C5: PD bif at 2nd PD (double), kernel=2 squares, compose stripes by hand 
aux.sw=-1;aux.coeff=[1 1];p=poswibra('2dh1','bpt2','2dpd2a',ds,aux); % stripes 
p.nc.dsmax=0.5; p.hopf.bisec=2; pause; p=cont(p,30); 
aux.sw=-1;aux.coeff=[1 0];p=poswibra('2dh1','bpt2','2dpd2b',ds,aux); % squares 
p.nc.dsmax=0.5; p.hopf.bisec=2; pause; p=cont(p,30); 
%% C6: BD of b-cont, || ||
fnr=3; figure(fnr); clf; c=8; plotbra('hom2Db','pt10',fnr,c,'cl','k','fp',6,'lp',12); 
plotbra('2dh1','pt47',fnr,c,'cl','b'); 
plotbra('2dpd1','pt45',fnr,c,'cl','r','lab',20); 
plotbra('2dpd2a','pt34',fnr,c,'cl',p2pc('b1'),'lab',20); 
plotbra('2dpd2b','pt34',fnr,c,'cl','m','lab',30); 
ylabel('max(|u_1,u_2|)'); 
axis([9.9 10.6 3.3 8.7]); 
%% C7: sol plots  
v=[0 90]; aux.lay=[1,4]; aux.pind=[1 15 29 43];
mhf('2dpd1','pt20',1,2,aux); title('2dpd1/pt20, u1(x0,t)'); xlabel('t'); pause 
mhf('2dpd2a','pt20',1,1,aux); title('2dpd2a/pt20, u1(x0,t)'); xlabel('t'); pause
mhf('2dpd2b','pt30',1,1,aux); title('2dpd2b/pt30, u1(x0,t)'); xlabel('t'); 