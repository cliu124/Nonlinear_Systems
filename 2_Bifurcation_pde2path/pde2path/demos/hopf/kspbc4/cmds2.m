%% 2ndary bif to traveling wave 
p=swibra('1','bpt1','1-1',-0.1); p.file.smod=10; p.nc.bisecmax=20; 
p.u0x=p.mat.M0*(p.mat.Kx*p.u(1:p.nu)); p=cont(p,80); 
%% secondary Hopf from traveling wave -> modulated TW 
aux.tl=50; aux.xif=0.1; ds=0.1; aux.pcfac=10; aux.dlam=0; aux.nqnew=0; 
p=hoswibra('1-1','hpt1',ds,4,'h4',aux); p=belon(p,4);
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=0; p.hopf.fltol=1e-4; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=0.3; p.nc.dsmin=1e-4; p.file.smod=2; p.hopf.y0dsw=1; 
p.nc.tol=1e-6; p=cont(p,20); 
%% Multipliers 
[muv1,muv2,ind,h]=floqpsap('h4','pt5'); pause; [muv1,muv2,ind,h]=floqap('h4','pt20'); 
%% plot BD, here speed over alpha 
pcmp=4; figure(3); clf; 
plotbra('1','pt75',3,pcmp,'cl',p2pc('b1'),'fp',60);
plotbra('1-1','pt30',3,pcmp,'cl',p2pc('r2'),'lp',30);
plotbra('h4','pt20',3,pcmp,'cl','g','lab',20);
xlabel('\alpha'); ylabel('s');
%% soln plots  
mclf(1); plotsol('1-1','hpt1'); title('1-1/hpt1'); pause 
hoplotf('h4','pt10',1,1); figure(1); title('h4/pt20'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); shading interp;
