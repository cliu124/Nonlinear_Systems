%% bru on (rather small 1D) domain 
p=[]; lx=2*pi/0.7; nx=50; a=3; b=8; du=6; dv=10; 
par=[a b du dv]; p=bruinit(p,lx,nx,par); 
p=initeig(p,10); p.nc.neig=[4, 4]; 
%% continue in D_u (to Turing) 
p.nc.ilam=3; p=setfn(p,'hom1Da'); p.sol.ds=-0.1; p=findbif(p,2); p=cont(p,10); 
%% stripes
p=swibra('hom1Da','bpt1','T1',0.1); p.nc.lammin=3; p=cont(p,20); 
%% continue T1 in b (stays stable) 
p=swiparf('T1','pt9','T1b',2); p.sol.ds=0.1; p.sw.verb=2; p.usrlam=11; p=cont(p,12); 
%% BD and soln plots of stripes 
fnr=3; figure(fnr); clf; c=8; plotbra('hom1Da','pt20',fnr,c,'cl','k'); 
plotbra('T1','pt15',fnr,c,'cl','k','lab',9); ylabel('||u_1||'); pause 
clf(fnr); plotbra('T1b','pt11',fnr,c,'cl','k','lab',11); ylabel('||u_1||'); pause 
plotsol('T1','pt9'); pause; plotsol('T1b','pt11');
%% analytic Hopf & Turing lines
du=2.6:0.2:7; bT=(1+a*sqrt(du./dv)).^2; bH=1+a^2+0*du; figure(1); clf; plot(du,bT,'-k'); hold on; 
plot(du,bH,'-b'); legend('T','H'); xlabel('D_u'); ylabel('b'); 
text(6.5,8,'1'); text(3,8,'2'); text(3,11,'3'); text(6.5,11,'4'); axis tight;
%% init on small (arbitrary) 1D domain, and use spufu to plot disp rel. 
p=[]; lx=2*pi/0.7; nx=50; a=3; b=8; du=6.5; dv=10; om=0.42; al=0; 
par=[a b du dv al om]; p=bruinit(p,lx,nx,par); p.nc.dsmax=0.5; 
p.fuha.ufu=@spufu; p0=p; p.sol.ds=1e-6; p.sol.dsmax=1e-6; p=setfn(p,'dummy'); p=cont(p,1);
p.u(p.nu+1:p.nu+4)=[3;8;3;dv]; p=cont(p,1);
p.u(p.nu+1:p.nu+4)=[3;11;3;dv]; p=cont(p,1); p.u(p.nu+1:p.nu+4)=[3;11;6.5;dv]; p=cont(p,1);
%% HP cont and BP cont; init at small b, d_u to cross Turing line in b, then cont Turing in d_u
% for HP cont use HP from hom1db 
p=[]; lx=2*pi/0.7; nx=50; a=3; b=6; du=3; dv=10; om=0.42; al=0; 
par=[a b du dv al om]; p=bruinit(p,lx,nx,par); p.nc.dsmax=0.5; p.sw.jac=1; 
p.nc.ilam=2; p.sol.ds=0.1; p=setfn(p,'dummy'); p=cont(p,5); 
%%
p=bpcontini('dummy','bpt1',3,'bpl'); huclean(p); 
p.sol.ds=0.1; p.nc.del=1e-2; p.plot.bpcmp=2;  
p.sw.spjm=1; p.sw.spjac=1; p.sw.jac=1; % spjac=0 OK, but jac=0 not working yet 
[Ja,Jn]=bpjaccheck(p); Jd=abs(abs(Ja-Jn)); e1=max(max(Jd)); mclf(10); spy(Jd>0.4*e1); pause 
p=cont(p,10); 
%%
p=hpcontini('hom1Db','hpt1',3,'hpc1'); p.nc.del=1e-4; p.sol.ds=-0.1;
[Ja,Jn]=hpjaccheck(p); % check hpjac
p.plot.bpcmp=2;  p=cont(p,10); 
figure(3); clf; plotbra('hpl','pt10',3,2,'cl','b'); 
plotbra('bpl','pt10',3,2,'cl','k'); box on; 
%%
hpcontexit('hpc1','pt10','du1'); pause; figure(2); clf; aux=[]; aux.tl=40; 
p=hoswibra('du1','hpt1',0.1,4,'duh1',aux); p.nc.dsmax=0.2; p=setbel(p,2,1e-4,5,@lss); 
p.plot.bpcmp=8; p=cont(p,10); 
%%
p=bpcontini('dummy','bpt1',3,'bpl'); huclean(p); p.sol.ds=0.1; p.nc.del=1e-2; 
p.plot.bpcmp=2; [Ja,Jn]=bpjaccheck(p); 
Jd=abs(Ja-Jn); e1=max(max(Jd)); mclf(6); spy(Jd>e1*0.9); 
%% 
p=hpcontini('hom1Db','hpt1',3,'hpl'); p.nc.del=1e-4; [Ja,Jn]=hpjaccheck(p); 
Jd=abs(Ja-Jn); e1=max(max(Jd)); mclf(6); spy(Jd>e1*0.9); 









