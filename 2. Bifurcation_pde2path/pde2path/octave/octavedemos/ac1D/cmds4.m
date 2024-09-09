%% Demo for imperfect pitchforks 
close all; keep pphome; 
%% c2: init (generic), then specific settings (could also be set in init)
p=[]; par=[1 -0.2 1 0.1 0]; p=acinit(p,5,30,par); 
ulam=[0 0.5]; p.usrlam=ulam; p=setfn(p,'i1'); cont(p,20);
%% init at larger lamda to find other branches 
p=[];  par=[1 0.2 1 0.1 0]; p=acinit(p,5,20,par);
p.usrlam=ulam; p=setfn(p,'i2a');p.sol.ds=0.05; p=cont(p,20);
p=[]; par=[1 0.2 1 0.1 0]; p=acinit(p,5,20,par); 
p.usrlam=ulam; p=setfn(p,'i2b'); p.sol.ds=-0.05; p=cont(p,20);
p=[]; par=[1 0.6 1 0.1 0]; p=acinit(p,5,20,par); 
p.usrlam=ulam; p=setfn(p,'i3a'); p=cont(p,20);
p=[]; par=[1 0.6 1 0.1 0]; p=acinit(p,5,20,par); 
p.usrlam=ulam; p.sol.ds=-0.05; p=setfn(p,'i3b'); p=cont(p,20);
%% branch plots
figure(3); clf; f=3; c=0; plotbra('i1',f,c,'cl','r'); 
plotbra('i2a',f,c,'cl','b'); plotbra('i2b',f,c,'cl','b'); 
plotbra('i3a',f,c,'cl','m'); plotbra('i3b',f,c,'cl','m'); 
%% solution plots 
plotsol('i1','pt4',1,1,1); pause; % using pause to export figure
plotsol('i1','pt10',1,1,1); pause; plotsol('i1','pt18',1,1,1); pause; 
plotsol('i2b','pt4',1,1,1); pause; plotsol('i2b','pt12',1,1,1); pause 
plotsol('i2a','pt15',1,1,1); pause; plotsol('i3b','pt0',1,1,1); pause 
plotsol('i3b','pt12',1,1,1); 
%% different solutions using deflation (here at fixed lambda=0.1) 
p=[]; par=[1 0 1 0.1 0]; p=acinit(p,5,100,par); p.nc.dsmax=0.05; 
[u,res,iter]=nloop(p,p.u); p.u=u; plotsol(p); % compute and store 1st solution 
x=getpte(p); x=x'; x2=max(x); % point coordinates (to set Iguesses for defl)
amps=[0.2 0.25 1 -1]; % amplitudes for Iguesses 
ug=p.u; global p2pglob; al3=1; p=deflinit(p,al3); % init deflation solver 
p.defl.nsw=2; % reset selected deflation parameters/switches, here norm 
p.nc.tol=1e-4; p.defl.al3=1; p.sw.newt=2; % newt=0,1: nloop, 2=NLEQ1, 3=fsolve
for i=1:length(amps) % deflation loop
    ug(1:p.nu)=amps(i)*cos(pi/(2*x2)*x); % set initial guess
    [u, p]=deflsol(p,ug); plotsol(p); pause  % run deflated Newton, plot sol
end 
p1=p; % store solutions from deflation 
%% run cont on solutions from deflation
p=postdefl(p1,1,'d1a',0.1); p=cont(p,20); % 1st branch, positive direction 
p=postdefl(p1,2,'d2a',0.1); p=cont(p,20); % 2nd branch, positive direction 
p=postdefl(p1,1,'d1b',-0.1); p=cont(p,20);
p=postdefl(p1,2,'d2b',-0.1); p=cont(p,20);
%% branch plots 
f=3; c=0; figure(f); clf;
plotbra('d1a',f,c,'cl','r'); plotbra('d1b',f,c,'cl','r');
plotbra('d2a',f,c,'cl','b'); plotbra('d2b',f,c,'cl','b');
axis([-0.2 0.51 0.1 3.5]); 
%% solution plots 
plotsol('d1a','pt10',1,1,1); pause; plotsol('d2a','pt10',1,1,1); 