%% commands for 1D KS, incl. Hopf branches 
close all; keep pphome; 
%% C1: init and zero-branch 
p=[]; lx=2; nx=100; % domain and discr
al=0.42; m=0; par=[al; m; 0; 0];  % m=mass, par(3)=eps, par(4)=s (speed) 
p=ksinit(p,nx,lx,par); p=setfn(p,'0'); screenlayout(p);
p.nc.mu1=10; p.nc.mu2=1; % large spacing of evals, be loose about localization  
p.nc.ilam=[1 3]; p=cont(p,40); % initial steps 
p.sol.ds=-0.001; p.nc.dsmax=0.001; % some more steps with smaller stepsize
p.nc.mu2=5; p=cont(p,20); p.file.smod=10; p=cont(p,10); 
%% C2: compute Turing-branches 
for i=1:3; 
    is=mat2str(i); p=swibra('0',['bpt' is],is,i*0.01); 
    p.sw.bifcheck=0; p=cont(p,5); % a few steps without PC 
    p.file.smod=20; p.u0x=p.mat.M0*(p.mat.Kx*p.u(1:p.nu)); % set reference profile for phase-invariance 
    p.nc.nq=2; p.nc.ilam=[1 3 4]; p.tau=[p.tau; 0]; % now switch on PCs 
    p.sw.bifcheck=2; p.nc.dsmax=0.2; p.nc.tol=1e-6; p=cont(p,i*70); 
end 
%% C3: 1st Hopf bifurcation 
figure(2); clf; ds=0.1; clear aux; aux.dlam=0; aux.nqnew=0; aux.tl=30; 
aux.xif=0.1; aux.y0dsw=0; p.hopf.nfloq=20; 
aux.nqh=2; aux.qfh=@qfh; aux.qfhder=@qfhjac; % func handles to hopf constraints 
p=hoswibra('2','hpt1',ds,4,'h1',aux); p.nc.ilam=1; p.hopf.ilam=[3 4]; 
p.hopf.fltol=5e-2; p.hopf.nfloq=10; p.hopf.flcheck=2; p.sw.verb=0; 
p.hopf.sec=1; p.nc.tol=1e-4; p.nc.dsmax=0.3; p.file.smod=2; p=cont(p,1); 
p.hopf.flcheck=1; p=cont(p,14); % floqps fails for larger amplitudes, 
% hence switch to floq: caution, only large multipliers seem correct
%% 2nd Hopf, with finer t-discret 
aux.tl=50; aux.xif=1; ds=0.1; p=hoswibra('3','hpt1',ds,4,'h2',aux); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=1; p.hopf.fltol=1e-1; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=0.3; p.nc.dsmin=1e-4; p.file.smod=2; p.hopf.y0dsw=2; 
p.nc.tol=1e-4; p=cont(p,30); 
%% 3rd Hopf 
aux.tl=50; aux.xif=0.5; ds=0.1; p=hoswibra('4','hpt1',ds,4,'h3',aux); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=1; p.hopf.fltol=1e-1; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=0.5; p.file.smod=2; p.hopf.y0dsw=2; 
p.nc.tol=1e-1; p=cont(p,2); % start with very large tol, gradually decrease 
p.nc.tol=1e-2; p=cont(p,2); p.nc.tol=1e-3; p=cont(p,6); 
p.nc.tol=1e-6; p=cont(p,20); % decrease tol once we're on the Hopf branch
%% Multipliers 
[muv1,muv2,ind,h]=floqpsap('h1','pt5'); pause; 
[muv1,muv2,ind,h]=floqpsap('h2','pt15'); pause
[muv1,muv2,ind,h]=floqpsap('h3','pt30');
%% 2ndary bif to traveling wave 
p=swibra('1','bpt1','1-1',0.1); p.file.smod=10; 
p.u0x=p.mat.M0*(p.mat.Kx*p.u(1:p.nu)); p=cont(p,40); 
%% global BD plot 
pcmp=6; figure(3); clf; plotbra('0','pt70',3,pcmp,'cl','k','lsw',0);
plotbra('1','pt75',3,pcmp,'cl',p2pc('b1'));
plotbra('2','pt100',3,pcmp,'cl',p2pc('b2'));
plotbra('3','pt100',3,pcmp,'cl',p2pc('b3')); 
plotbra('1-1','pt30',3,pcmp,'cl',p2pc('r1'),'lab',30,'fp',1);
xlabel('\alpha'); ylabel('max(u)');
%% local BD plot 
pcmp=6; figure(3); clf; %pcmp=3; 
plotbra('2','pt100',3,pcmp,'cl',p2pc('b2')); %,'lab',30,'fp',1);
plotbra('3','pt160',3,pcmp,'cl',p2pc('b3')); %,'lab',50,'fp',1);
plotbra('4','pt160',3,pcmp,'cl',p2pc('b4'));% ','lab',50,'fp',1);
plotbra('h1','pt15',3,pcmp,'cl','r','lab',[5,15]);
plotbra('h2','pt30',3,pcmp,'cl','m','lab',30);
plotbra('h3','pt30',3,pcmp,'cl',p2pc('br1'),'lab',30);
axis([0.005 0.06 9 25]); xlabel('\alpha'); ylabel('max(u)'); 
%% sol plots 
plotsol('2','hpt1'); title('2/hpt1'); xlabel([]); pause; plotsol('3','hpt1'); xlabel([]); title('3/hpt1'); pause
plotsol('4','hpt1'); title('4/hpt1'); pause
%%
hoplotf('h1','pt5',1,1); figure(1); title('h1/pt5'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause;
hoplotf('h1','pt15',1,1); figure(1); title('h1/pt15'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause;
hoplotf('h2','pt30',1,1); figure(1); title('h2/pt30'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause
hoplotf('h3','pt30',1,1); figure(1); title('h3/pt30'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); 
%% secondary Hopf from traveling 
aux.tl=100; aux.xif=0.1; ds=0.1; aux.pcfac=10; clf(2); 
p=hoswibra('1-1','hpt1',ds,4,'h4',aux); p.nc.imax=10; 
%p.u0x=p.mat.Kx*p.u(1:p.nu); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.flcheck=0; p.hopf.fltol=1e-1; p.sw.verb=2;
p.hopf.sec=1; p.nc.dsmax=0.3; p.nc.dsmin=1e-4; p.file.smod=2; p.hopf.y0dsw=2; 
p.nc.tol=1e-3; p=cont(p,5); % start with very large tol, gradually decrease 
p.nc.tol=1e-4; p=cont(p,15); 
%% Multipliers 
[muv1,muv2,ind,h]=floqpsap('h4','pt10'); pause; [muv1,muv2,ind,h]=floqpsap('h4','pt20'); 
%% plot BD, here speed over alpha 
pcmp=4; figure(3); clf; 
plotbra('1','pt75',3,pcmp,'cl',p2pc('b1'),'fp',60);
plotbra('1-1','pt40',3,pcmp,'cl',p2pc('r2'),'lp',30);
plotbra('h4','pt20',3,pcmp,'cl','g','labi',10);
xlabel('\alpha'); ylabel('s');
%% soln plot 
plotsol('1-1','hpt1'); title('1-1/hpt1'); pause
hoplotf('h4','pt10',1,1); figure(1); title('h4/pt10'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1');  shading interp; pause;
hoplotf('h4','pt20',1,1); figure(1); title('h4/pt20'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); shading interp;
%% check tint
p=loadp('h1','pt5'); p.u(1:p.nu)=p.hopf.y(1:p.nu,10); 
t0=0; ts=[]; dt=0.005; nt=50; nc=0; pmod=10; smod=0; plotsol(p); 
%%
[p,t0,ts,nc]=tintxs(p,t0,ts,dt,nt,nc,pmod); 