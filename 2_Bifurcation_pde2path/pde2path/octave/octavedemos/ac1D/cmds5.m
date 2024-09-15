%% Demo for imperfect pitchforks 
close all; keep pphome; 
%% c2: init (generic), then specific settings (could also be set in init)
p=[]; par=[1 -0.2 1 0.1 0]; p=acinit(p,5,40,par); p.nc.dsmax=0.05; 
ulam=[0 0.5]; p.usrlam=ulam; p=setfn(p,'i1b'); cont(p,10);
%% use deflation
p=loadp('i1','pt6'); x=getpte(p); x=x'; x2=max(x); % plotsol(p); pause;
p.osol=p.u; u0(1:p.nu)=p.u(1:p.nu)+0.2*cos(pi/(2*x2)*x); % set initial guess for Newton 
% p.u=u0; plotsol(p); pause 
p.sw.jac=1; p.deflq=2; p.deflal1=1e-3; p.deflal2=1; %p.nc.del=1e-6; 
p.nd=1; [u, p]=deflsol(p,u0); p.u=u; plotsol(p); pause  % run deflated Newton, use result for next deflation 
u0(1:p.nu)=p.u(1:p.nu)+1*cos(pi/(2*x2)*x); [u, p]=deflsol(p,u0); p.u=u; % pause 
u0(1:p.nu)=p.u(1:p.nu)-1.5*cos(pi/(2*x2)*x); [u, p]=deflsol(p,u0); % 4th solution 
p.nc.dsmax=0.1; p1=p; % save p in p1 to start cont. runs
%% run cont on obtained solutions
p=postdefl(p1,2,'j2a',0.1); p=cont(p,20);
p=postdefl(p1,2,'j2b',-0.1); p=cont(p,20);
p=postdefl(p1,4,'j3a',0.1); p.sw.bifcheck=0; p=cont(p,20);
p=postdefl(p1,4,'j3b',-0.1); p.sw.bifcheck=0; p=cont(p,20);
%% 
f=3; c=0; figure(f); clf;
plotbra('j2a',f,c,'cl','r'); 
plotbra('j2b',f,c,'cl','r');
plotbra('j3a',f,c,'cl','b'); 
plotbra('j3b',f,c,'cl','b');
axis([-0.2 0.51 0.1 3.5]); 
%% jaccheck 
p=loadp('i1','pt6'); x=getpte(p); x=x'; x2=max(x); % plotsol(p); pause;
p.osol=p.u; u0(1:p.nu)=p.u(1:p.nu)+0.2*cos(pi/(2*x2)*x)+0.05; p.u=u0; 
p.nc.del=1e-6; p.deflq=2; p.deflal1=1e-2; p.deflal2=1; p.nd=1; 
p.fuha.sGb=p.fuha.sG; p.fuha.sGjacb=p.fuha.sGjac;  % mod rhs
p.fuha.sG=@deflsG; p.fuha.sGjac=@deflsGjac; 
[Gu, Gun]=jaccheck(p); Gud=abs(Gu-Gun); e1=max(max(Gud)); 
spy(Gud>e1/1000); 
p.fuha.sG=p.fuha.sGb; p.fuha.sGjac=p.fuha.sGjacb; 


