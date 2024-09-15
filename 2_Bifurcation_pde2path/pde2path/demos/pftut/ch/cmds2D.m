close all; keep pphome; p=[]; 
%% CH 2D, eps=0.1, init and cont of homogen. branch 
p=[]; m=-1.2; eps=1/10; lam=0; par=[m eps lam]; lx=[0.5 0.5]; nx=[30 30]; 
p=chinit(p,lx,nx,par); p=setfn(p,'2D'); p.nc.dsmax=0.04; p=cont(p,100); 
%% 1st BP, double, hence cswibra 
p0=cswibra('2D','bpt1',aux); p0.nc.lammax=0.3; 
p=seltau(p0,1,'2D1-sp',3); p.sol.ds=0.1; p=cont(p,130); % spots
%% use gentau for stripes 
p=gentau(p0,[1 0.5],'2D1-st',3); p.sol.ds=-0.1; p=cont(p,60); % stripes
%% 2ndary bifs 
p=swibra('2D1-sp','bpt2','2D-2a',-0.1); p=cont(p,10); 
%% 2nd BP, simple 
p=swibra('2D','bpt2','2D2',0.1); p=cont(p,90); 
%% 3rd BP, double again 
aux=[]; p0=cswibra('2D','bpt3',aux); 
p=seltau(p0,2,'2D3-sp',3); p.sol.ds=0.1; p=cont(p,50); 
p=gentau(p0,[0.4 0.4],'2D3-st'); p.sol.ds=0.1; p=cont(p,40); 
%% BD plot E over m
figure(3); clf; plotbra('2D','pt60','lsw',0); 
plotbra('2D1-sp','pt80','cl','b','lab',[30 60]); 
plotbra('2D1-st','pt60','cl','r','lab',40); 
plotbra('2D2','pt90','cl','m','lab',50); 
plotbra('2D3-sp','pt50','cl',p2pc('o1'),'lab',20); 
plotbra('2D3-st','pt40','cl',p2pc('o2'),'lab',10); 
axis([-0.9 0.1 0.5 2.8]); ylabel('E_\epsilon'); 
%% BD plot E over m, zoom into 1-spots near m=0; 
figure(3); clf; plotbra('2D1-sp','pt80','cl','b','lab',[55 60 65 70]); 
plotbra('2D-2a','pt10','cl',p2pc('g1'),'lab',5); 
ylabel('E_\epsilon'); axis([-0.1 0.1 1.17 1.27]); 
%% sol plots
mypsol('2D1-sp','pt30'); mypsol('2D1-sp','pt55'); mypsol('2D1-st','pt40'); 
mypsol('2D2','pt50'); mypsol('2D3-sp','pt20');mypsol('2D3-st','pt10'); 
mypsol('2D1-sp','pt60'); mypsol('2D1-sp','pt65'); mypsol('2D-2a','pt5');
%% continuation in eps of FP1 on 2D1-sp: 
figure(2); clf; p=spcontini('2D1-st','fpt1',2,'fp2c'); p.sol.ds=-0.001; plotsol(p); 
p.sw.bifcheck=0; p.nc.dsmax=0.05; p.file.smod=5; p.usrlam=[0.075 0.05]; 
p.nc.del=1e-4; p.nc.tol=1e-6; huclean(p); p=cont(p,20); 
%% BD for fold-cont 
figure(3); clf; plotbra('fp2c','pt20',3,1,'lab',[0 10]); xlabel('\epsilon'); 
%% soln plots 
mypsol('fp2c','pt0'); mypsol('fp2c','pt20');
%% FP cont exit 
huclean(p); p=spcontexit('fp2c','pt12','sp2a'); p.sw.spcalc=1; p.usrlam=0; 
p.nc.dsmax=0.2; p.nc.lammax=1; p.nc.lammin=-1; p.plot.bpcmp=5;  p=cont(p,30); 
%% BD at small eps 
figure(3); clf; plotbra('sp2a','pt40','cl','b'); 
axis([-0.4 0.4 1.1 1.4]); ylabel('E_\epsilon'); 
%% trulle MA (just illustration) 1741-> 1151 
p=loadp('sp2a','pt26'); plotsol(p,10,1,1);view(20,40); nola; p.np, pause 
p.sw.trul=1; op=troptions2D(); op.etafac=5e-6; %u=@etafua;
op.verbose=2; op.innerit=2; 
p.trop=op; p.trop.sw=15; p=oomeshada(p,'ngen',5); plotsol(p,10,1,1); nola; view(20,40); 
%% trulle MA,  1741-> 1559
p=loadp('sp2a','pt68', 'sp2r'); plotsol(p); p.np, 
p.sw.trul=1; op=troptions2D(); op.etafac=5e-6; op.verbose=2; op.innerit=2; 
p.trop=op; p.trop.sw=15; p=oomeshada(p,'ngen',5); plotsol(p,10,1,1); nola; view(20,40); 
%%
[E1, E22]=chE(p,p.u); E2,E22
%% cont with trulle MA 
p.nc.amod=10; p=cont(p,100); 
%% plot BD 
figure(3); clf; plotbra('sp2','pt100','cl','b');
%% cont in eps to see if E_eps->|I|
p=swiparf('2D1-st','pt40','2D1-st-eps',[2 3]); getaux(p)', p.nc.lammax=1; clf(2); 
p.sol.ds=-0.01; p.nc.dsmax=0.04; p.nc.tol=1e-8; p=cont(p,10); 