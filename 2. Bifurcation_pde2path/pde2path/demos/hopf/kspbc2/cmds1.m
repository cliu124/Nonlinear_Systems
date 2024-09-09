%% commands for 1D KS, incl. Hopf branches 
close all; keep pphome; 
%% init and zero-branch 
p=[]; lx=2; nx=70; % domain and discr
al=0.42; m=0; par=[al; m; 0; 0];  % m=mass, par(3)=eps, par(4)=s (speed) 
p=ksinit(p,nx,lx,par); p=setfn(p,'0'); screenlayout(p);
p.nc.mu1=10; p.nc.mu2=5; % large spacing of evals, be loose about localization  
p.sol.ds=-0.001; p.nc.dsmax=0.001; % some more steps with smaller stepsize
p.nc.ilam=[1 3]; p.file.smod=100; p.nc.nsteps=1e4; p=findbif(p,3); 
%%
p=cont(p,100); % initial steps 
p.sol.ds=-0.001; p.nc.dsmax=0.001; % some more steps with smaller stepsize
p.nc.mu2=5; p=cont(p,20); p.file.smod=10; p=cont(p,10); 
%% Turing-branches 
for i=1:3; 
    is=mat2str(i); p=swibra('0',['bpt' is],is,i*0.01); 
    p.file.smod=50; p.sw.bifcheck=0; p.nc.tol=1e-6; p=cont(p,5); % a few steps without PC 
    %u0x=p.mat.M(1:p.nu/2,1:p.nu/2)*(p.mat.Kx*p.u(1:p.nu/2));
    u0x=p.mat.Kx*p.u(1:p.nu/2); 
    p.u0x=u0x;  % set profile for transl-invariance 
    p.nc.nq=2; p.nc.ilam=[1 3 4]; p.tau=[p.tau; 0]; % now switch on PCs 
    p.sw.bifcheck=2; p.nc.bisecmax=20; p.sw.verb=0; p.plot.pmod=20; 
    p.nc.dsmax=i/2; p.nc.tol=1e-6; p=cont(p,round(i^1.8)*100); 
end 
%% 1st Hopf bifurcation  
figure(2); clf; ds=2; clear aux; aux.dlam=0; aux.nqnew=0; aux.bpcmp=6; aux.tl=40; 
aux.xif=0.1; %aux.y0dsw=1;  %aux.pcfac=1000; 
aux.nqh=2; aux.qfh=@qfh; aux.qfhder=@qfhjac; % func handles to hopf constraints 
p=hoswibra('2','hpt1',ds,4,'h1',aux); p.nc.ilam=1; p.hopf.ilam=[3 4]; 
p.hopf.fltol=1e-2; p.hopf.nfloq=10; p.hopf.flcheck=1; p.sw.verb=0; 
p.hopf.sec=0; p.nc.tol=1e-6; p.nc.dsmax=2; p.file.smod=2; 
p=belon(p,3); % Hopf-borders are size 4: PC, arclength, qh (twice) 
p=cont(p,20); % floqps fails for larger amplitudes, 
%% 2nd Hopf, with finer t-discret 
aux.tl=50;  ds=0.5; aux.y0dsw=1; aux.pcfac=100; 
p=hoswibra('3','hpt1',ds,4,'h2',aux); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.fltol=1e-2; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=5; p.nc.dsmin=1e-4; p.file.smod=2; p=belon(p,3);
p.nc.tol=1e-4; p=cont(p,10); p.nc.tol=1e-6; p=cont(p,20); 
%% 3rd Hopf 
aux.tl=60; aux.xif=0.1; ds=2; aux.pcfac=100; 
p=hoswibra('4','hpt1',ds,4,'h3',aux); 
p.nc.ilam=1; p.hopf.ilam=[3 4]; p.hopf.fltol=1e-2; p.sw.verb=0;
p.hopf.sec=1; p.nc.dsmax=10; p.file.smod=2;  p=belon(p,3);
p.nc.tol=1e-4; p=cont(p,10);  % start with large tol, gradually decrease 
p.nc.tol=1e-6; p=cont(p,10); 
%% Multipliers 
[muv1,muv2,ind,h]=floqap('h1','pt10'); pause; 
[muv1,muv2,ind,h]=floqpsap('h2','pt20'); pause
[muv1,muv2,ind,h]=floqap('h3','pt20');
%% global BD plot 
pcmp=6; figure(3); clf; plotbra('0','pt70',3,pcmp,'cl','k','lsw',0);
plotbra('1','pt100',3,pcmp,'cl',p2pc('b1'));
plotbra('2','pt300',3,pcmp,'cl',p2pc('b2'));
plotbra('3','pt600',3,pcmp,'cl',p2pc('b3')); 
%plotbra('1-1','pt30',3,pcmp,'cl',p2pc('r1'),'lab',30,'fp',1);
xlabel('\alpha'); ylabel('max(u)');
%% local BD plot 
pcmp=6; figure(3); clf; %pcmp=3; 
plotbra('2','pt300',3,pcmp,'cl',p2pc('b2')); %,'lab',30,'fp',1);
plotbra('3','pt600',3,pcmp,'cl',p2pc('b3')); %,'lab',50,'fp',1);
plotbra('4','pt800',3,pcmp,'cl',p2pc('b4'));% ','lab',50,'fp',1);
plotbra('h1','pt20',3,pcmp,'cl','r','lab',20);
plotbra('h2','pt30',3,pcmp,'cl','m','lab',20);
plotbra('h3','pt20',3,pcmp,'cl',p2pc('br1'),'lab',20);
%axis([0.005 0.06 9 25]); xlabel('\alpha'); ylabel('max(u)'); 
%% sol plots 
plotsol('2','hpt1'); title('2/hpt1'); xlabel([]); pause; plotsol('3','hpt1'); xlabel([]); title('3/hpt1'); pause
plotsol('4','hpt1'); title('4/hpt1'); pause
%%
hoplotf('h1','pt10',1,1); figure(1); title('h1/pt10'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause;
hoplotf('h2','pt20',1,1); figure(1); title('h2/pt20'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); pause
hoplotf('h3','pt20',1,1); figure(1); title('h3/pt20'); view(20,70); 
xlabel('x'); ylabel('t'); zlabel('u_1'); 