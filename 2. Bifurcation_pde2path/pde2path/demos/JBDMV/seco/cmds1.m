%% MWBS97 model, trivial and Turing branches, and bif from Turing 
al=0.02; j0=3.4; tau=0.05; D=8; kc=(al*tau/D)^0.25; nw=8; % parameters, 
lx=nw*pi/kc; nx=nw*40; par=[j0;al;tau;D]; dir='0'; % dom.size and nx
p=secoinit(lx,nx,par); p=setfn(p,dir); % init, set output dir for tr.branch
%% compute guess for Hopf-spectral shift, then cont trivial branch 
p=initeig(p,0.5); p.nc.neig=[10,10]; p.nc.eigref(1)=-0.05; 
p.nc.dsmax=0.05; p.sw.verb=1; p.file.smod=5;  p=cont(p,30); 
%% primary Turing mode 
p=swibra('0','bpt1','T1'); p.nc.dsmax=0.05; p.sol.ds=-0.02; p=cont(p,60);
%% T-mode via cswibra, useful for eigenvector inspection 
aux.m=6; aux.besw=0; p0=cswibra('0','bpt2',aux);
%% now generate bif directions and go 
p=gentau(p0,[1],'T2'); p.nc.dsmax=0.1; p=cont(p,40);
p=gentau(p0,[0 0 1],'T3'); p.nc.dsmax=0.1; p=cont(p,40);
%% T4-mode via cswibra
aux.m=6; aux.besw=0; p0=cswibra('0','bpt3',aux); pause 
p=gentau(p0,-1,'T4'); p.nc.dsmax=0.1; p=cont(p,40);
%% T-front snake 
clf(1); p=swibra('T1','bpt1','T1-1'); p.nc.dsmax=0.05; p.sw.foldcheck=1; 
p.sol.ds=-0.02; p.file.smod=10; p.nc.mu2=0.01; p.sw.jac=0; 
tic; p=cont(p,100);toc 
%% HP10 on front-snake 
para=4; ds=0.1; aux=[]; aux.tl=30; aux.dlam=0; 
p=hoswibra('T1-1','hpt10',ds,para,'T1-1-H10',aux); pause  
p.hopf.flcheck=0; p.sw.verb=0; 
p=setbel(p,2,1e-4,5,@lss); p.nc.dsmax=0.5; p=cont(p,50); 
%% HP18 
p=hoswibra('T1-1','hpt18',ds,para,'T1-1-H18',aux); pause  
p.hopf.flcheck=0; p.sw.verb=0; 
p=setbel(p,2,1e-4,5,@lss); p.nc.dsmax=0.5; p=cont(p,50); 
%% loc.Turing  snake 
clf(1); p=swibra('T1','bpt2','T1-2'); p.nc.dsmax=0.05; p.sw.foldcheck=1; 
p.file.smod=10; p.sol.ds=-0.02; p=cont(p,200);
%% HP8 on loc.Turing snake 
para=4; ds=0.1; dsmax=1.1; xi=1e-3; aux.tl=20; 
p=hoswibra('T1-2','hpt8',ds,para,'T1-2-H8',aux); pause  
p.hopf.flcheck=0; p.sw.verb=0; p=setbel(p,2,1e-4,5,@lss); p.nc.dsmax=0.5; p.plot.bpcmp=8; 
p.file.smod=10; p=cont(p,100); 
%% BD-plot,
fnr=3; figure(fnr); clf; cmp=8; plotbra('0',fnr,cmp,'cl','k','lsw',0); 
plotbra('T1','pt60',fnr,cmp,'cl','b','lab',20); 
plotbra('T2','pt40',fnr,cmp,'cl',p2pc('b2'),'lsw',0); 
plotbra('T3','pt40',fnr,cmp,'cl',p2pc('g1'),'lsw',0); 
plotbra('T4','pt40',fnr,cmp,'cl',p2pc('g2'),'lab',25); 
plotbra('T1-1','pt400',fnr,cmp,'cl',p2pc('o1'),'lab',350,'lsw',0,'lp',380,'hplab',[10 18]); 
plotbra('T1-1-H10','pt100',fnr,cmp,'cl',p2pc('r1'),'lsw',0,'lab',[20 50]);
plotbra('T1-1-H18','pt100',fnr,cmp,'cl',p2pc('r2'),'lsw',0); 
axis([3.05 3.7 0 1.8]); xlabel('j_0'); ylabel('||u_1||_*'); 
%% BD plot, zoom of loc T-snake 
figure(fnr);clf; 
plotbra('T1-2','pt100',fnr,cmp,'cl',p2pc('o2'),'lsw',0,'hplab',8,'lp',175); 
plotbra('T1-2-H8','pt100',fnr,cmp,'cl',p2pc('r2'),'lsw',0,'lab',[20 80],'lp',90);
axis([3.24 3.35 1.15 1.45]); xlabel('j_0'); ylabel('||u_1||_*'); 
%% soln plots, steady 
figure(1); clf; 
plotsol('T1','pt30',1,1,1); noltix; pause; 
plotsol('T4','pt35',1,1,1); noltix; pause; 
plotsol('T1-1','hpt10',1,1,1); noltix; pause; 
plotsol('T1-1','hpt14',1,1,1); noltix; pause; 
plotsol('T1-1','hpt18',1,1,1); noltix; pause; 
plotsol('T1-1','pt370',1,1,1); noltix
%% hoplots, with Fl-multipl., uncomment the desired directory 
v=[15 60]; aux.nfloq=20; 
dir='T1-1-H10'; ptli=[20 50]; 
dir='T1-2-H8'; ptli=[20 80]; 
for i=ptli; %120 % 10:10:130;
    pt=['pt' mat2str(i)]; hoplotf(dir,pt,1,1); figure(1); 
    title([dir '/' pt]); view(v); shading interp; ylabel('t'); 
    muv1=floqap(dir,pt,aux);  pause
end
