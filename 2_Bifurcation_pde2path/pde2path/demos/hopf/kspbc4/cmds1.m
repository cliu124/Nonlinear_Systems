%% commands for 1D KS, incl. Hopf branches 
close all; keep pphome; 
%% init and zero-branch 
al=0.42; m=0; par=[al; m; 0; 0];  % m=mass, par(3)=eps, par(4)=s (speed) 
p=[]; lx=2; nx=200;p=ksinit(p,nx,lx,par); p=setfn(p,'0'); % domain and discr
p.nc.mu1=10; p.nc.mu2=1; % large spacing of evals, be loose about localization  
p.nc.ilam=[1 3]; p=cont(p,40); % initial steps 
p.sol.ds=-0.001; p.nc.dsmax=0.001; % some more steps with smaller stepsize
p.nc.mu2=5; p=cont(p,20); p.file.smod=10; p=cont(p,10); 
%% compute branches of steady patterns
for i=1:4
    is=mat2str(i); p=swibra('0',['bpt' is],is,i*0.01); 
    p.file.smod=20; p.sw.bifcheck=0; p=cont(p,5); % a few steps without PC 
    p.u0x=p.mat.Kx*p.u(1:p.nu); % set profile for transl-invariance 
    p.nc.nq=2; p.nc.ilam=[1 3 4]; p.tau=[p.tau; 0]; % now switch on PCs 
    p.sw.bifcheck=2; p.nc.dsmax=0.2; p.nc.tol=1e-6; 
    p.plot.pmod=20; 
    p=cont(p,30+i*40); 
end 
%% 1st Hopf bifurcation 
figure(2); clf; ds=0.1; clear aux; aux.dlam=0; aux.nqnew=0; aux.tl=30; aux.bpcmp=6;
aux.xif=0.1; aux.y0dsw=1; % use PDE to set d/dt u_0 for phase-constr. (in t)
aux.nqh=2; aux.qfh=@qfh; aux.qfhder=@qfhjac; % func handles to hopf constraints 
p=hoswibra('2','hpt1',ds,4,'h1',aux); p.nc.ilam=1; p.hopf.ilam=[3 4]; 
p.hopf.fltol=1e-2; p.hopf.nfloq=10; p.hopf.flcheck=1; p.sw.verb=0; 
p=belon(p,3); % Hopf-borders are size 4: PC, arclength, qh (twice) 
p.hopf.sec=1; p.nc.tol=1e-6; p.nc.dsmax=0.3; p.file.smod=5; p=cont(p,1);
p.hopf.flcheck=1; p=cont(p,14); % floqps fails for larger amplitudes, 
% hence switch to floq: caution, only large multipliers seem correct
%% uniform mesh-ref in t: increase # time-slices by fac, then cont 
%p=loadp('h1','pt5','h1ref'); fac=1.5; p=uhopftref(p,fac); p=cont(p,5); 
%% 2nd Hopf, with finer t-discret 
aux.tl=40; aux.xif=10; ds=0.1; aux.pcfac=100; 
p=hoswibra('3','hpt1',ds,4,'h2',aux); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=0; p.hopf.fltol=1e-2; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=0.5; p.nc.dsmin=1e-4; p.file.smod=2; %p.hopf.y0dsw=2; 
p=belon(p,3); 
p.nc.tol=1e-4; p=cont(p,10); % first with large tol 
p.nc.tol=1e-6; p=cont(p,10); % decrease tol 
%% 3rd Hopf 
aux.tl=40; aux.xif=0.5; ds=0.1; aux.pcfac=100; 
p=hoswibra('4','hpt1',ds,4,'h3',aux); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=1; p.hopf.fltol=1e-1; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=0.5; p.file.smod=2; p.hopf.y0dsw=1; 
p=belon(p,3); 
p.nc.tol=1e-2; p=cont(p,2); p.nc.tol=1e-3; p=cont(p,6); 
p.nc.tol=1e-6; p=cont(p,24); % decrease tol once we're on the Hopf branch
%% Multipliers 
[muv1,muv2,ind,h]=floqap('h1','pt10'); pause; 
[muv1,muv2,ind,h]=floqpsap('h2','pt4'); pause
[muv1,muv2,ind,h]=floqap('h3','pt30');
%% global BD plot 
pcmp=6; figure(3); clf; plotbra('0','pt70',3,pcmp,'cl','k','lsw',0);
plotbra('1','pt80',3,pcmp,'cl',p2pc('b1'));
plotbra('2','pt100',3,pcmp,'cl',p2pc('b2'));
plotbra('3','pt200',3,pcmp,'cl',p2pc('b3')); 
plotbra('1-1','pt80',3,pcmp,'cl',p2pc('r1')); 
xlabel('\alpha'); ylabel('max(u)');
%% local BD plot 
pcmp=6; figure(3); clf; pcmp=6; plotbra('2','pt100',3,pcmp,'cl',p2pc('b2')); 
plotbra('3','pt150',3,pcmp,'cl',p2pc('b3')); plotbra('4','pt150',3,pcmp,'cl',p2pc('b4')); 
plotbra('h1','pt15',3,pcmp,'cl','r','lab',[5,15]);
plotbra('h2','pt20',3,pcmp,'cl','m','lab',50,'lp',50);
plotbra('h3','pt30',3,pcmp,'cl',p2pc('br1'),'lab',30);
axis([0.005 0.06 9 25]); xlabel('\alpha'); ylabel('max(u)'); 
%% steady sol plots 
plotsol('2','hpt1'); title('2/hpt1'); xlabel([]); pause; 
plotsol('3','hpt1'); xlabel([]); title('3/hpt1'); pause
plotsol('4','hpt1'); title('4/hpt1'); 
%% Hopf soln plots
hoplotf('h1','pt5',1,1); figure(1); title('h1/pt5'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); shading interp; pause;
hoplotf('h1','pt15',1,1); figure(1); title('h1/pt15'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1');  shading interp; pause;
hoplotf('h2','pt20',1,1); figure(1); title('h2/pt20'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1');  shading interp; pause
hoplotf('h3','pt30',1,1); figure(1); title('h3/pt30'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); shading interp; 