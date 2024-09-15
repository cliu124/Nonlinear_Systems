%% demo breather, init 
close all; format compact; keep pphome; 
p=[]; al=0.11; bet=1; ga=6; del=0.01; D=2; s=0; par=[al bet ga del D s]; 
lx=50; pw=lx/4; nx=80; p=brinit(p,lx,nx,par,pw); plotsol(p); 
x=getpte(p);u=1./cosh(0.1*x');v=u;p.u0=[u;v];p.u0x=p.Dx*p.u0; % reference profile
%% converge to TW (here standing pulse) via freezing
nt=1e4; pmod=250; vmod=100; dt=0.25; vel=[];  t1=0; 
[p,t1,vel]=tintfreeze(p,t1,dt,nt,pmod,vel,vmod); 
%% gradually decrease del=par(4) to first find breather at delta approx 0.006, 
% which becomes unstable near del=0.0044
% for cont-setting, go directly to next cell
p.u(p.nu+4)=p.u(p.nu+4)-0.0005; p.u(p.nu+4), t1=0; vel=[]; 
[p,t1,vel]=tintfreeze(p,t1,dt,nt,pmod,vel,vmod); 
figure(7); clf; h=plot(vel(1,:),vel(2,:)); axis tight; p0=p; 
%% set up and continue extended stat system; 
p=setfn(p,'s1'); p.u0=p.u(1:p.nu); p.u0x=p.Dx*p.u0; % reset reference profile 
p.nc.nq=1;p.sw.qjac=1; p.fuha.qf=@qf;p.fuha.qfder=@qfder; % switch on constraints
p.nc.ilam=[4 6]; % set active pars to delta and s 
figure(2); clf; p.sol.ds=-1e-3; p.sw.spcalc=1; p=cont(p,20); % continue 
%% swibra to Hopf 
aux.nqnew=0; aux.dlam=0; % switch off steady constr., use trivial init. predictor
aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac; % switch on average speed constr. 
aux.tw=1e-3; % small weight of T in arclength useful here 
%aux.tl=30; indir='s1'; dir='h1b'; aux.xif=1000; fltol=5e-2; % quick 
aux.tl=70; indir='s1'; dir='h1'; aux.xif=1000; fltol=5e-2; % more accurate  
% use rather large fltol to ignore spurious (tl too small) bif with ga=1
p=hoswibra(indir,'hpt1',0.25,4,dir,aux); p.hopf.flcheck=1; p.hopf.jac=1; 
p=belon(p,2); % switch on bel (bordered elimination) 
p.u0=p.u(1:p.nu); p.u0x=p.Dx*p.u0; % reset reference profile to prof.at HBP 
p.nc.ilam=4; p.hopf.ilam=6; p.file.smod=1; p.nc.dsmax=2; p.hopf.fltol=fltol; 
p.nc.tol=1e-6; p.sw.verb=0; p.hopf.bisec=5; p=cont(p,30); 
%% more bisections yield better localization of PD point (but bisec=5 also works) 
% this in particular needed for tl>50 
p=loadp('h1','pt19'); p.sol.ds=0.1;  p.nc.dsmax=0.1; 
p.hopf.fltol=5e-2; p.hopf.bisec=5; p=cont(p,10); 
%% period doubling from breather
huclean(p); ds=0.1; p=poswibra('h1','bpt1','pd1',ds); p.nc.tol=0.01; p.nc.dsmax=2; 
p=belon(p,1); 
p.hopf.nqh=0; p.hopf.ilam=[]; % switch off average speed constraints (allows better stepsizes) 
p.sw.bifcheck=0; p.hopf.flcheck=1; p.nc.dsmin=1e-3; p.sw.verb=0; p.hopf.fltol=1e-2; 
p=cont(p,2); p.nc.tol=1e-2; p=cont(p,5); p.nc.tol=1e-6; p=cont(p,23); 
%% branch plotting
bf=3; pcmp=8; figure(bf); clf; plotbra('s1','pt20',bf,pcmp,'cl','k','fp',2, 'lp',12);
plotbra('h1','pt30',bf,pcmp,'cl','r','lab',[15,30]); 
xlabel('\delta'); ylabel('max(u_1)'); 
%%
bf=3; pcmp=8; figure(bf); clf; 
plotbra('h1','pt30',bf,pcmp,'cl','r','fp',12); 
plotbra('pd1','pt30',bf,pcmp,'cl','m','lab',[10 20],'fp',1); 
xlabel('\delta'); ylabel('max(u_1)');
%% soln plots
dir='h1'; pt='pt'; ind=[15 30]; 
for i=1:length(ind);
    q=loadp(dir, [pt mat2str(ind(i))]); 
    hoplot(q, 1,1,1); figure(1); view(10,70); xlabel('x'); ylabel('t'); 
    title(['u_1 at ' dir '/' [pt mat2str(ind(i))]]); shading interp;
    set(gca,'ZTick',[0 0.4 0.8]); figure(6); axis tight; pause 
end 
dir='pd1'; pt='pt'; ind=[10 20]; 
for i=1:length(ind);
    q=loadp(dir, [pt mat2str(ind(i))]); 
    hoplot(q, 1,1,1); figure(1); view(10,70); xlabel('x'); ylabel('t'); 
    title(['u_1 at ' dir '/' [pt mat2str(ind(i))]]); shading interp;
    set(gca,'ZTick',[0 0.4 0.8]); figure(6); axis tight; pause 
end 
%% Floquet checks and plots 
[muv1,muv2,ind,mo,h]=floqap('h1','pt15'); pause
[muv1,muv2,ind,mo,h]=floqap('h1','pt30'); pause
[muv1,muv2,ind,mo,h]=floqap('pd1','pt10'); pause 
[muv1,muv2,ind,mo,h]=floqap('pd1','pt20'); 