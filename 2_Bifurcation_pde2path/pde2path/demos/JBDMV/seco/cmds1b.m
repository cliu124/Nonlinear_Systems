%% BP and HP cont, needs a larger delta for FDs for the param-derivatives
p=bpcontini('0','bpt1',2,'bpc1'); p.sol.ds=0.005; p.nc.dsmax=0.2; 
p.nc.del=1e-3; p.plot.bpcmp=1; p.nc.tol=1e-6; p=cont(p,50); 
%% reverse direction 
p=loadp('bpc1','pt5','bpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% primary Hopf
p=hpcontini('0','hpt1',2,'hpc1');  p.sol.ds=0.05; p.nc.dsmax=0.2;  
p.nc.del=1e-2;  p.fuha.spjac=@hpjac; %[Ga,Gn]=hpjaccheck(p); 
p.plot.bpcmp=1; p.nc.tol=1e-3; p=cont(p,20); 
p=loadp('hpc1','pt5','hpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% standard branch plot 
mclf(4); plotbra('bpc1',4,1,'cl','b'); plotbra('hpc1',4,1,'cl','r'); 
axis([0 0.1 1 3.5]); ylabel('j0'); xlabel('\alpha'); 
%% plotbradat 
mclf(3); p=loadp('bpc1'); h1=plotbradat(p,3,7,8,'-b'); hold on; 
set(h1,'linewidth',2); p=loadp('hpc1'); h3=plotbradat(p,3,7,8,'d-r');
p=loadp('bpc1b'); h1=plotbradat(p,3,7,8,'-b'); set(h1,'linewidth',2); 
p=loadp('hpc1b'); h3=plotbradat(p,3,7,8,'d-r');
axis([1 3.6 0 0.1]); legend('T_1','H_1'); xlabel('j0'); ylabel('\alpha'); 
annotation('textarrow',[0.64 0.58],[0.62 0.62],'String','(c)','Fontsize',14);
annotation('textarrow',[0.84 0.84],[0.4 0.28],'String','(d)','Fontsize',14);
text(3.22,0.05,'(b)','Fontsize',14); set(gca,'fontsize',14); 
%% MWBS, disp.relations, uncomment the desired par-values in the next line 
al=0.04; j0=3.3; %al=0.06; j0=2.4; %al=0.01; j0=3.5; 
tau=0.05; D=8; kc=(al*tau/D)^0.25; nw=8; % parameters, 
lx=nw*pi/kc; nx=nw*40; par=[j0;al;tau;D]; dir='0b'; % dom.size and nx
p=secoinit(lx,nx,par); p=setfn(p,dir); % init, set output dir for tr.branch
p.fuha.ufu=@spufu; p.sol.ds=1e-6; p=cont(p,1); figure(10); 
title(['(j0,\alpha)=(' mat2str(j0) ',' mat2str(al) ')']); 