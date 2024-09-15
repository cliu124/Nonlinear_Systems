close all; keep pphome; p=[]; 
%% CH 3D, init and cont of hom 
m=-0.6; eps=sqrt(1/20); lam=0; par=[m eps lam]; lx=[0.5 0.5 0.5]; n=25; nx=[n n n]; 
p=chinit(p,lx,nx,par); p=setfn(p,'3D'); p.nc.nq=1; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1; p.nc.ilam=[1,3]; % aux eqns 
p.nc.dsmax=0.04; p.nc.mu2=0.05; p.sw.verb=2; p.nc.tol=1e-6; p=cont(p,25); 
%% 1st BP, triple, hence cswibra, also set usrlam for continuation
aux=[]; aux.m=3; aux.besw=1; p0=cswibra('3D','bpt1',aux); p0.nc.lammax=0.05; 
p0.pm.resfac=1e-3; p0.usrlam=[-0.3 0 0.3]; p0.sol.ds=0.02; p0.sw.bifcheck=0; 
%%
p=seltau(p0,11,'3Dl',3); p.nc.dsmax=0.1;  p=cont(p,50); 
p=seltau(p0,7,'3Dtu',3); p=pmcont(p,30); % the tubes need pmcont 
p=seltau(p0,14,'3Dsp',3); p.nc.dsmax=0.02; p=cont(p,15); 
p.nc.dsmax=0.1; p.nc.lammax=0.05; p=cont(p,20); 
%% BD plot E over m
figure(3); clf; plotbra('3D','pt25','lsw',0); 
plotbra('3Dsp','pt30','cl','b','lab',[15 27],'fp',1); 
plotbra('3Dl','pt17','cl','r','lab',15,'fp',1); 
plotbra('3Dtu','pt30','cl','m','lab',[17],'fp',1);
axis([-0.5 0.01 0.7 1.2]); ylabel('E_\epsilon'); 
%% sol plots
v=[7,30]; mypsol('3Dl','pt15',v); mypsol('3Dtu','pt17',v); 
plotsol('3Dsp','pt15',1,1,2); nolti; zticks([]); zlabel(''); axis 'image'; view(v); pause; 
plotsol('3Dspr','pt27',1,1,2); nolti; zticks([]); zlabel(''); axis 'image'; view(v); 
%% cont in eps, sp15 
p=swiparf('3Dsp','pt15','3Dsp-eps',[2 3]); p.sol.ds=-0.01; clf(2); p.nc.lammax=1; 
p.nc.dsmax=0.1; p=cont(p,20); 
%% BD and sol plot 
figure(3); clf; plotbra('3Dsp-eps','pt10',3,5,'lab',10); 
xlabel('\epsilon'); ylabel('E_\epsilon'); 
plotsol('3Dsp-eps','pt10',1,1,2); nolti; zticks([]); zlabel(''); axis 'image'; view(v);
%% cont in eps, sp at m=0
p=swiparf('3Dsp','pt27','3Dsp0-eps',[2 3]); getaux(p)', p.nc.lammax=1; 
p.sol.ds=-0.01; p.nc.dsmax=0.04; p.nc.tol=1e-8; p=cont(p,20); 
%% cont in eps, lamella at m=0,  
p=swiparf('3Dl','pt15','3Dl-eps',[2 3]); getaux(p)', p.nc.lammax=1; 
p.sol.ds=-0.01; p.nc.dsmax=0.04; p.nc.tol=1e-8; p=cont(p,20); 
%% BD plot 
figure(3); clf; plotbra('3Dsp0-eps','pt20',3,5,'lab',15,'lp',15); 
plotbra('3Dl-epsr','pt15',3,5,'cl','r','lab',15); 
xlabel('\epsilon'); ylabel('E_\epsilon');
%% sol-plot
plotsol('3Dsp0-eps','pt15',1,1,2,'sh',0); nolti; zticks([]); zlabel(''); axis 'image'; view(v);
%% cont in eps, lamella at m=0, with trulle
p=swiparf('3Dl','pt15','3Dl-epsr',[2 3]); getaux(p)', p.nc.lammax=1; pause 
p.sol.ds=-0.01; p.nc.dsmax=0.04; p.nc.tol=1e-8; 
p.sw.trul=1; op=troptions3D(); op.etafac=5e-5; op.verbose=2; op.innerit=2; 
p.trop=op; p.trop.sw=15; p=oomeshada(p,'ngen',3);  stansavefu(p); 
%%
p.nc.amod=5; p=cont(p,20); 
%% sol-plot
plotsol('3Dl-epsr','pt15',1,1,3,'sh',0); nolti; zticks([]); zlabel(''); axis 'image'; view(v);
%%
[E1,E2]=chE(p,p.u); 
%%
p=loadp('3Dl-epsr','pt15'); p.np