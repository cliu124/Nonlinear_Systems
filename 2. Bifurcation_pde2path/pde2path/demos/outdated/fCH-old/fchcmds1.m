%% fCH, commands for straight (initial) interfaces 
%% initial stripe 
close all; keep pphome; p=[];
lx=pi/2; ly=2*pi; p.pp=3; p.del=0; sw=1; wi=0.1; a=1.25; 
nx=60; ny=120; p.fuha.wfu=@wfu2; % potential and derivatives 
p.eta1=1; p.eta2=3; p.eps=0.35; p.mi=-0.9; % p.mi=-0.9;
p=fchinits(p,lx,ly,nx,ny,wi,a,sw); p.nc.tol=1e-10;p.nc.lammax=2; 
%% try to converge! repeat this cell to pick suitable mesh by hand!
p.nc.ngen=1; p=meshref(p,'maxt',20000);p=setbmesh(p);
%% now continue 
p=pmcont(p);
%% Bif from first 4 primary bif., and one 2ndary 
p=swibra('p','bpt1','q1',0.01); p=pmcont(p); 
p=swibra('p','bpt2','q2',0.01); p=pmcont(p);
p=swibra('p','bpt3','q3',0.01); p=pmcont(p);
%% plot BD
figure(3);clf;cmp=3;
plotbra('p',3,cmp,'cl','k','ms',0,'lwun',3,'lwst',3);
plotbra('q1',3,cmp,'cl','r','lab',30,'lwun',3,'lwst',3);
plotbra('q2',3,cmp,'cl','m','lwun',3,'lwst',3);
plotbra('q3',3,cmp,'cl','b','lab',30,'lwun',3,'lwst',3);
xlabel('\eta_1'); ylabel('\gamma');title(['mass=-0.9']);
axis([1 2 -20.2 -18.6]);
%% solution plots 
plotsol('p','bpt3',1,1,2); xlabel(''); ylabel(''); pause 
plotsol('q1','pt20',1,1,2);  xlabel(''); ylabel(''); 
%% tangent plot via calling swibra again, then some settings
q=swibra('p','bpt3','du',0.01);xlabel(''); ylabel(''); colorbar off;
title('\tau_1 at p/bpt3','FontSize',p.plot.fs);
set(gca,'FontSize',q.plot.fs); 
%% pearling from iguess / use sfem=0 for this!, 
lx=pi/2; ly=3*pi/4; p.pp=3; p.del=0; sw=2; wi=0.1; a=1.25; 
nx=50; ny=50; p.fuha.wfu=@wfu2; % potential and derivatives; nx=60; ny=60;
p.eta1=1; p.eta2=3; p.eps=0.35; p.mi=-0.98; 
p=fchinits(p,lx,ly,nx,ny,wi,a,sw); p0=p; p=setfn(p,'pp'); p0=p;
%% again, repeat this cell as appropriate
% p.fsol.fsol=1; % fsolve a good alternative for this 
p.sw.sfem=0; p.nc.ngen=1; p=meshref(p,'maxt',15000); p=setbmesh(p);
%% cont works well on this one 
p=cont(p); 

