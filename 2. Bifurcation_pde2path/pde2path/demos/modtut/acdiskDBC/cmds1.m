%% ac2D disk, DBCs, chebychev-setup 
close all; keep pphome; 
%% init 
p=[]; s=0; par=[0.5 -0.1 0 -1 s 0]; % c,lam,quad,cubic; coeff of bdry term 
lx=5; ly=1; nr=16; nphi=32; p=acinit(p,lx,nr,nphi,par); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=0.5; p.sw.verb=2; p.nc.neig=20; p.nc.lammax=1; p.nc.eigref=-1;  
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8; p.plot.bpcmp=7;  p=cont(p,10);
%% radial branch 
p=swibra('tr','bpt1','b1',0.1); p=cont(p,10); 
%% angular branches, need PC 
ind=[2:4]; 
for i=ind; 
  is=mat2str(i); p=swibra('tr',['bpt' is],['b' is],0.1); p=cont(p,5); pause 
  p.file.smod=10; p.nc.ilam=[2; 5];  p.nc.nq=1;  
  p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; % analytical jac for aux. eqn.
  p.sw.bprint=[2 5]; p.nc.dsmax=0.11;  p=cont(p,15); 
end  
%% BD plot 
f=3; c=7; figure(f); clf;
plotbra('tr','pt10',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab',6); 
plotbra('b2',f,c,'cl','r','lab',6); 
plotbra('b3',f,c,'cl','m','lsw',0); 
plotbra('b4',f,c,'cl',p2pc('g1')); 
ylabel('||u||_2'); axis([0 1.01 0 0.8]); 
%% soln plots 
plotsol('b1','pt6'); pause; plotsol('b2','pt6'); pause; 
plotsol('b4','pt3'); 
%%
p=loadp('b2','pt6'); p.ups=1; p.ipz=1; p.uview=[20,70]; plotsol(p); 
%%
p=loadp('b4','pt3'); p.ups=4; p.ipz=1;  plotsol(p); 