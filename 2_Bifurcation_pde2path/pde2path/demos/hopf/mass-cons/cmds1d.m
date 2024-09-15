%% mass-cons, demo for POs in mass-conserving system, 
% see 'symtut/modfro' and other symtut-demos for 'hopf-constraints'  
close all; format compact; keep pphome; 
%% C1: init, and continuation of hom branch 
ndim=1; dir='hom1d'; p=[]; lx=pi; nx=100; % domain size and spat.resolution 
par=[1; 10; 0; 0; 0; 1];  % d1 d2 d3 m a1 b1
p=mcinit(p,lx,nx,par,ndim); p=setfn(p,dir); % initialization 
p.nc.nq=1; p.fuha.qf=@qf;  % 1 steady constraint (mass), and its func.handle
p.sw.qjac=1; p.fuha.qfder=@qfjac; % use analytical jac for q, and func.handle 
p.nc.xiq=0.1; p.nc.ilam=[5 6]; % weight of constraint in arclength, active vars 
p.nc.dsmax=0.2; p=cont(p,50);  % run the continuation 
%% C2: hopf with constraints, passed to hoswibra via aux vars in aux 
para=4; ds=0.005; aux=[]; aux.dlam=0; aux.nqnew=0; aux.tl=50; 
aux.xif=100; aux.pcfac=100;  % weight factors, see hostanparam  
aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac; % func.handles to hopf contraints 
for i=2:2 % continue for 15 steps, first three with large tol
  p=hoswibra('hom1d',['hpt' mat2str(i)],ds,para,['h' mat2str(i)],aux); 
  p.nc.ilam=5; p.hopf.ilam=6; p.sw.verb=0; p.hopf.sec=1; p.nc.dsmax=0.5; 
  p.file.smod=5; p.nc.tol=1e-4; p=belon(p,2); p.sw.bifcheck=0; p.hopf.flcheck=0; 
  p=cont(p,3); p.nc.tol=1e-6; p=cont(p,22); 
end
%% C3: time integrate from some point on Hopf orbit, preparations 
p=loadp('h2','pt15'); p.u(1:p.nu)=p.hopf.y(1:p.nu,1); % load Hopf point, and 
% adjust some settings; in particular store system stiffness matrix in p.mat.K 
dir='t1'; u0=p.u(1:p.nu); p=setfn(p,dir); p.plot.shading='interp'; 
ts=[]; t0=0; npp=50; nt=200; pmod=50; smod=2; tsmod=1; nc=0; 
par=p.u(p.nu+1:end); Ks=p.mat.K; p.mat.K=[[par(1)*Ks par(2)*Ks];[par(3)*Ks Ks]]; 
%% C4: time-integration (repeat if necessary) 
[p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); p0=p; 
%% C5: x-t plot (of both components)  
si=0; incr=4; vv=[30,70]; nt=50; 
tintplot1d(dir,si,incr,nt,1,1,vv); shading interp; xlabel('x'); ylabel('t'); title('u1(x,t)');
tintplot1d(dir,si,incr,nt,2,2,vv); shading interp; xlabel('x'); ylabel('t'); title('u2(x,t)');
%% C6: continue from result of tint; again hopf for decreasing alpha 
p=p0; plotsol(p); p.sol.restart=1; p.sol.ds=-0.1; p.sw.para=0; % reset settings to 
% steady case, in particular restore scalar stiffness matrix (-Laplacian)
p.mat.K=Ks; p=rmfield(p,'hopf'); p=setfn(p,'hom1d2'); p=resetc(p); 
p.nc.nq=1; p.nc.ilam=[5 6]; p.fuha.headfu=@stanheadfu; p.fuha.ufu=@stanufu; 
p=cont(p,10); 
%% plot BD, 1-6=pars, 7=T, 8,9=max/min(u1), 10=L^2(u1), 11=mass, 12,13: max/min(u2)
bpcmp=9; wnr=3; figure(wnr); clf
switch bpcmp
    case 6; yl='\beta'; case 9;yl='min(u_1)'; case 11; yl='mass'; case 12;yl='max(u_2)'; 
end 
plotbra('hom1d','pt50',3,bpcmp,'lsw',3,'fp',6,'lp',50); % label only HBPs 
plotbra('h1','pt25',3,bpcmp,'lab',15,'cl','b');   
plotbra('h2','pt25',3,bpcmp,'lab',15,'cl','r');   
plotbra('h3','pt25',3,bpcmp,'lab',15,'cl','m');  
xlabel('\alpha'); ylabel(yl); 
%% plot solns 
v=[15,50]; 
hoplotf('h1','pt15',1,1); figure(1); title('u1 at h1/pt15'); shading interp; 
xlabel('x'); ylabel('t'); zlabel('u_1'); view(v); pause;
hoplotf('h1','pt15',1,2); figure(1); title('u2 at h1/pt15'); shading interp; 
xlabel('x'); ylabel('t'); zlabel('u_2'); view(v); pause;
hoplotf('h2','pt20',1,1); figure(1); title('u1 at h2/pt20'); shading interp; 
xlabel('x'); ylabel('t'); zlabel('u_1'); view(v); pause; 
hoplotf('h3','pt15',1,1); figure(1); title('u1 at h3/pt15'); shading interp; 
xlabel('x'); ylabel('t'); zlabel('u_1'); view(v); pause
hoplotf('h4','pt10',1,1); figure(1); title('h4/pt15'); shading interp; 
xlabel('x'); ylabel('t'); zlabel('u_1'); view(v);
%% plot Floquets
[muv1,muv2,ind,mo,h]=floqap('h1','pt15');