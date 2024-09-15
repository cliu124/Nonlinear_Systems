%% commands for computing Eck(k) for 1D SH 
keep pphome; close all; p=[]; 
%% init and zero-branch 
lx=20*pi; nx=round(20*lx); ly=0.1; ny=1; ndim=1; lam=-0.02; nu=0; l=1; 
par=[lam; nu; l]; p=shinit(p,nx,lx,ly,ndim,par); p.file.smod=5; 
p=setfn(p,'0'); p.usrlam=0; p.nc.bisecmax=10; p.nc.dsmax=0.02; p=findbif(p,3);
%% Turing-branches 1-3 (wave-numbers 1, 19/20 and 21/20 (for lx=20*pi))
% resp. 19/20 and 21/20 for lx=20*pi
%huclean(p); 
p=swibra('0','bpt1','b1',0.01);  p=cont(p,20); 
p=swibra('0','bpt2','b2',0.01);  p=cont(p,20);
p=swibra('0','bpt3','b3',0.01);  p=cont(p,20);
%% 2ndary branches
p=swibra('b2','bpt1','2-1'); pause; p=cont(p,80); 
p=swibra('b3','bpt1','3-1'); pause; p=cont(p,20); 
p=swibra('b3','bpt3','3-2'); pause; p=cont(p,20); 
%% BD plot 
figure(3); clf; cmp=5; wnr=3; 
plotbra('b1','pt20',wnr,cmp,'cl','k','lsw',0); plotbra('b2','pt20',wnr,cmp,'cl','r'); 
plotbra('b3','pt20',wnr,cmp,'cl','m'); plotbra('2-1','pt20',wnr,cmp,'cl','r');
plotbra('3-1','pt20',wnr,cmp,'cl','m'); plotbra('3-2','pt20',wnr,cmp,'cl','m');
axis([0 0.03 0 0.1]); xlabel('\lambda'); ylabel('||u||'); box on; 
%%
plotsol('3-1','pt10'); axis tight; pause; 
plotsol('b2','pt10'); axis tight; pause; plotsol('b3','pt10'); axis tight; pause 
%% continuation of BP1 from 0 in l (to obtain 'existence of stripes' region) 
figure(2); clf; p=bpcontini('0','bpt1',3,'bpc0a'); p.k=1; p.sol.ds=0.01; 
p.nc.lammin=-1; p.plot.bpcmp=1; p.sw.spjac=1; p.nc.del=1e-2; 
p.fuha.spjac=@bpjac; p.sw.spcalc=0; p.sw.bifcheck=0;
p.nc.tol=1e-8; p.file.smod=5; huclean(p); p0=p; 
%[Ja, Jn]=bpjaccheck(p); Jd=abs(Ja-Jn); e1=max(max(Jd)); mclf(10); spy(Jd>0.1*e1); 
p=cont(p,10); 
%%
p=p0; p.sol.ds=-p.sol.ds; p=setfn(p,'bpc0b'); p=cont(p,10); 
%% continuation in l of BP1 on b2: b2 gains stab. here; k0 is 19/20, i.e., 
% dilated wave, hence continue to kap<1 for 'further dilation' 
figure(2); clf; p=bpcontini('b2','bpt1',3,'bpc2a'); p.k=19/20; p.sol.ds=-0.005; 
p.plot.bpcmp=1; p.sw.spcalc=0; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.tol=1e-6; 
p.sw.spjac=1; p.fuha.spjac=@bpjac; huclean(p); p.sw.jac=1; p.sw.spjac=2; % for testing 
[Ja, Jn]=bpjaccheck(p); Jd=abs(Ja-Jn); e1=max(max(Jd)); mclf(10); spy(Jd>0.1*e1); pause 
p=cont(p,20); 
%% continuation in l of BP3 on b3: b3 gains stab. here; k0 is 21/20, i.e., 
% compressed wave, hence continue to kap>1 for 'further compression' 
figure(2); clf; p=bpcontini('b3','bpt3',3,'bpc3a'); p.k=21/20; p.sol.ds=0.01; 
p.plot.bpcmp=1; p.sw.spjac=1; p.sw.bifcheck=0; p.nc.dsmax=0.05; p.nc.tol=1e-6; 
p.fuha.spjac=@bpjac; p.sw.spcalc=0; huclean(p); p=cont(p,20); 
%% plot 'existence line' (BP continuation of prim.bif) and stability lines 
% (continuation of BPs where sideband branches become stable) 
% k=l*k_branch at pos.5 (due to param-doubling during BP cont, lam at pos 1) 
figure(3); clf; hold on; p=loadp('bpc0a','pt10'); plotbradat(p,3,5,1); 
p=loadp('bpc0b','pt10'); plotbradat(p,3,5,1); 
p=loadp('bpc2a','pt20'); aux.ps='-r'; plotbradat(p,3,5,1,aux); 
p=loadp('bpc3a','pt20'); aux.ps='-m'; plotbradat(p,3,5,1,aux); 
k=linspace(0.5,1.5,40); kap=k.^2-1;
figure(3); hold on; plot(k,kap.^2,'*k'); % analytical existence line
plot(k,3*kap.^2-kap.^3,'*r'); % analytical stability line
axis([0.5 1.5 0 1]); xlabel('k'); ylabel('\lambda'); set(gca,'fontsize',16); 
%% jaccheck for BC, 
[Ja, Jn]=bpjaccheck(p); Jd=abs(Ja-Jn); me=max(max(Jd)), figure(10); spy(Jd>me/100); 
%% BP cont exit
huclean(p); p=bpcontexit('bpc1','pt5','1D2b'); %p=cont(p,10); 
