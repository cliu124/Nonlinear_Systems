%% 2ndary bif to traveling wave 
p=swibra('1','bpt1','1-1',0.1); p.file.smod=10; p.nc.bisecmax=20; p.nc.tol=1e-4; 
u0x=p.mat.Kx*p.u(1:p.nu/2);
p=cont(p,80); 
%% secondary Hopf from traveling wave -> modulated TW 
aux.dlam=0; aux.nqnew=0; aux.tl=50; ds=-0.2; aux.pcfac=100; 
aux.xif=0.1; aux.y0dsw=2; % use PDE to set d/dt u_0 for phase-constr. (in t)
aux.nqh=2; aux.qfh=@qfh; aux.qfhder=@qfhjac; % func handles to hopf constraints 
p=hoswibra('1-1','hpt1',ds,4,'h4',aux); pause 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=0; p.hopf.fltol=1e-4; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=1; p.nc.dsmin=1e-4; p.file.smod=2; p.hopf.y0dsw=2; 
beliluon; 
p.nc.tol=1e-4; p=cont(p,15); % start with large tol, gradually decrease 
p.nc.tol=1e-6; p=cont(p,15); 
%% Multipliers 
[muv1,muv2,ind,h]=floqpsap('h4','pt10'); pause; [muv1,muv2,ind,h]=floqap('h4','pt20'); 
%% plot BD, here speed over alpha 
pcmp=4; figure(3); clf; 
plotbra('1','pt200',3,pcmp,'cl',p2pc('b1'),'fp',60);
plotbra('1-1','pt80',3,pcmp,'cl',p2pc('r2'));
plotbra('h4','pt30',3,pcmp,'cl','g','lab',[10 20]);
xlabel('\alpha'); ylabel('s');
%% soln plots  
plotsol('1-1','hpt1'); title('1-1/hpt1'); pause
hoplotf('h4','pt10',1,1); figure(1); title('h4/pt10'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1');  shading interp; pause;
hoplotf('h4','pt20',1,1); figure(1); title('h4/pt20'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); shading interp;