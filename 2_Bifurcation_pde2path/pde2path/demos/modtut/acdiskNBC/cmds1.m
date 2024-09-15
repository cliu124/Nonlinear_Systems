%% ac2D disk, NBCs, chebychev-setup 
close all; keep pphome; 
%% init 
p=[]; s=0; par=[0.5 -0.1 0 -1 s]; % c,lam,quad,cubic
lx=5; ly=1; nr=20; nphi=40; p=acinit(p,lx,nr,nphi,par); p=setfn(p,'tr'); p.plot.pstyle=-1;  
p.nc.lammax=20; p.nc.dsmax=0.5; p.sw.verb=2; p.nc.neig=20; p.nc.lammax=1; 
p.sw.bifcheck=2; p.nc.mu1=1; p.nc.tol=1e-8; p.plot.bpcmp=6;  p=cont(p,10);
%%
ind=[1,4]; 
for i=ind
   is=mat2str(i); p=swibra('tr',['bpt' is],['b' is],-0.1); p.sw.verb=2; pause; 
   p.sw.bifcheck=0; p=cont(p,10); 
end
%%
ind=[2 3]; 
for i=ind; 
  is=mat2str(i); p=swibra('tr',['bpt' is],['b' is],0.1); p.sw.bifcheck=0; p=cont(p,3); 
  p.file.smod=10; p.nc.ilam=[2; 5];  p.nc.nq=1;  
  p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; % analytical jac for aux. eqn.
  p.sw.bprint=[2 5]; p.nc.dsmax=0.11;  p=cont(p,15); 
end  
%% BD plot 
f=3; c=6; figure(f); clf;
plotbra('tr','pt10',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lsw',0); 
plotbra('b2',f,c,'cl','r', 'lab',12); 
plotbra('b3',f,c,'cl','m','lab',14); 
plotbra('b4',f,c,'cl',p2pc('g1'),'lsw',0); 
ylabel('||u||_2'); axis([0 1 0 1]); 
%% soln plots 
plotsol('b2','pt12'); pause; plotsol('b3','pt14'); pause; plotsol('b4','pt5'); 
%%
p=loadp('b2','pt12'); p.ups=1; p.ipz=1; plotsol(p); 
%%
p=loadp('b4','pt5'); p.ups=4; p.ipz=1;  plotsol(p); 
%%
figure(10); clf; spy(p.mat.L); 