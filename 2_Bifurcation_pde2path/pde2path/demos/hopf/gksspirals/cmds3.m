%% cont of 1-armed spiral as rel.equilibr., and bif to meandering spirals
% alpha=1, delta=0.25 (convenient). 
% increased spatial resolution to have accurate bifs to mRW (and tip-paths)
close all; format compact; keep pphome; 
%% note: nx=45; mRWs are slow!  
ndim=2; dir='tr-b2'; p=[]; lx=pi; nx=45; r=-0.25; al=1; del=0.25; s=0; 
par=[r al del s]; p=rotinit(p,nx,par); p=setfn(p,dir); p.nc.neig=50; 
p.sol.ds=0.01; p.nc.dsmax=0.01; p.nc.ilam=1; p.nc.bisecmax=10; p=cont(p,20); 
%% compute RWs via twswibra
spar=4; kwnr=1; % speed (here rot) parameter index, and w-nr 
aux=[]; p=twswibra('tr-b2','hpt2',spar,kwnr,'rw1c',aux); p.L=2*pi; 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; 
p.u0x=p.mat.M\(p.mat.Krot*p.u0); p.u(1:p.nu)=p.u(1:p.nu)+0.1*p.tau(1:p.nu); 
plotsolu(p,p.u0,6,1,3); plotsolu(p,5*p.u0x,11,1,3); p.sw.bifcheck=2; 
p.nc.nq=1; p.nc.ilam=[1;spar];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; p.sw.verb=2; p.usrlam=[1 2 3]; 
p.sw.bprint=4; clf(2); p.nc.dsmax=1; p.sol.ds=0.02; pause; 
p.nc.tol=1e-4; p=cont(p,2); p.nc.tol=1e-8; p=cont(p,20); 
%% profiles from cont in r 
pstyle=2; pt='pt15'; plotsol('rw1c',pt,1,1,pstyle); clplot(1); plotsol('rw1c',pt,2,2,pstyle); 
clplot(2); pause; l1=0.1; l2=0.1; levplot('rw1c','pt15',1,0,[l1 l2]); 
%% mRW; 
aux=[]; aux.dlam=0; aux.nqh=1; aux.nqnew=0; aux.tl=20; 
aux.qfh=@qfh; aux.qfhder=@qfhjac; 
p=hoswibra('rw1c','hpt2',0.2,4,'mrw1',aux); p.hopf.ilam=4; p.nc.ilam=1; 
bw=2; beltol=1e-4; belimax=5; droptol=1e-3; AMGmaxit=200; 
AMG=0; % set AMG=1 if ilupack available 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setilup(p,droptol,AMGmaxit); end 
p.sw.verb=0; p.nc.dsmax=0.21; p.file.smod=1; p.sw.bifcheck=0; p.hopf.flcheck=0; p.sw.verb=2; pause; 
p.hopf.y0dsw=0; p.hopf.pcfac=1e-3; % small pcfac is more robust; (Krot*u is poor) 
huclean(p); p.usrlam=[]; p.nc.tol=1e-4; p=cont(p,20); 
%% mRW; tl=40, needs ilupack 
aux.tl=40; p=hoswibra('rw1c','hpt2',0.2,4,'mrw1b',aux); p.hopf.ilam=4; p.nc.ilam=1; 
droptol=1e-2; AMGmaxit=500;  p=setilup(p,droptol,AMGmaxit); 
p.sw.verb=0; p.nc.dsmax=0.21; p.file.smod=1; p.sw.bifcheck=0; p.hopf.flcheck=0; p.sw.verb=2; pause; 
p.hopf.y0dsw=0; p.hopf.pcfac=1e-3; % small pcfac is more robust; (Krot*u is poor) 
huclean(p); p.usrlam=[]; p.nc.tol=1e-3; p=cont(p,20); 
%% BD of cont in r 
bpcmp=5; figure(3); clf; plotbra('rw1c','pt20',3,bpcmp,'lab',[15],'cl','m','fp',1,'lp',17); 
plotbra('mrw1','pt20',3,bpcmp,'lab',[5 25],'cl','b','fp',0); xlabel('r'); ylabel('T'); 
%% Flowers etc 
hoaux=[]; hoaux.lay=[1 3]; hoaux.pind=[1 3 5]; hoaux.view=[0 90];
hoaux.xtics=''; hoaux.ytics=''; aux.hoaux=hoaux; 
aux.mr=1; aux.pertol=0.025; aux.tol=1e-3; l1=0.1; l2=l1; dir='mrw1'; 
for i=5:20:25; 
    i, pt=['pt' mat2str(i,2)]; plottip(dir,pt,11,[l1 l2],aux); 
    hoplotft(dir,pt,1,1,hoaux); lfplottip(dir,pt,10,[l1 l2],aux); pause
end
%% period/freq plots
iv=1:1:20; liv=length(iv); T1=zeros(1,liv); s1=T1; T2=T1; om1=T1; om2=T1; rv=T1; 
for i=1:liv
    pt=['pt' mat2str(iv(i))];  p=loadp(dir,pt); 
    rv(i)=p.u(p.nu+1); T2(i)=p.hopf.T; om2(i)=1/T2(i); 
    s1(i)=p.u(p.nu+p.spar); T1(i)=2*pi/s1(i); om1(i)=1/T1(i); 
end
figure(1); clf; plot(rv,T1./T2,'-*'); legend('T_1/T_2'); 
set(gca,'FontSize',p.plot.fs); xlabel('r'); 
%% continue RW in diff constant (domain size) 
p=swiparf('rw1c','pt15','rw1d',[3 4]); clf(2); p.nc.neig=20; p.usrlam=[0.02 0.05 0.1 0.2]; 
p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.lammin=0.02; p.nc.neig=50; p=cont(p,30); 
%% alternative way of cont with reset of PC every 5 steps 
for i=1:4; p=cont(p,5); p.u0(1:p.nu)=p.u(1:p.nu); p.u0x=p.mat.M\(p.mat.Krot*p.u0); end 
%% BD of cont in delta
bpcmp=5; figure(3); clf; plotbra('rw1d','pt20',3,bpcmp,'lab',20); 
xlabel('\delta'); ylabel('T'); set(gca,'YTick',[2.07 2.075]); 
%%
plotsol('rw1d',pt,1,1,2); clplot(1); % profile from cont in delta
%% spectral plot 
pt='pt20'; p=loadp('rw1d',pt); wnr=10; [muv,ev]=plotspec(p,wnr); 2*pi/p.u(p.nu+p.spar)
% plotsolu(p,real(ev(:,1)),1,1,2); plotsolu(p,imag(ev(:,1)),2,1,2); % Efus 
%% zoom near imag axis 
figure(4); clf; plot(real(muv), imag(muv),'*'); axis([0 2 -0.4 10]); set(gca,'FontSize',16); 
%% spectral plot, at rw1c/hpt2
p=loadp('rw1c','hpt2'); wnr=10; [muv,ev]=plotspec(p,wnr); 
% plotsolu(p,real(ev(:,1)),1,1,2); plotsolu(p,imag(ev(:,1)),2,1,2); 
%% Floquet a posteriori, choose points of interest! 
aux.nfloq=20; [muv1,~,~,~,h]=floqap('mrw1b','pt25',aux); 