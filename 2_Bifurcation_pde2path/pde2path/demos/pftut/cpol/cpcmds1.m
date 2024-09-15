%%  cell polarization, Cussedu-Edelstein-Keshet-..-Madzvamuse 2018, eps=0.1, medium mehs 
close all; keep pphome; 
%% init 
p=[]; e=0.1; k0=0.1; ga=1; R=1; bpstyle=4; m=11; s=0; mlam=0; cut=[-2 0.01 -2]; 
par=[R  k0  e  ga  s  m  mlam]; % s=speed for x-PC, m=mass, mlam=lag-mult for m-cons
%    1                6
lx=pi; del=5e-2; ly=pi/2-del; 
nx=13; ny=10; ref=1; nr=8; h0=0.1; odir='h0'; 
p=cpinit(p,lx,ly,nx,ny,par,ref,bpstyle,cut,h0); 
p=setfn(p,odir); userplot(p,21); p.sw.para=2; 
%% checks of numjac, initial call to get sparsity pattern S
global pjS; pjS=0; 
r=resi(p,p.u); p.sw.jac=0; tic; Gu=getGu(p,p.u,r); toc, mclf(1), spy(Gu); title('sparsity of G_u'); 
tic; Gu2=getGu(p,p.u,r); toc; mclf(12); spy(Gu2); % subsequent call to see speedup 
%% but we run with jac=1 anyway; homogeneous branch 
p.nc.tol=1e-6; p.nc.ilam=[6 7]; p.nc.nq=1; p.sw.jac=1; 
p=cont(p,20);
%% primary pattern
p=swibra(odir,'bpt1','c1',0.01);  p.nc.tol=1e-6; p.sw.bifcheck=0; % switch off bifdetec for speed 
p.nc.dsmax=0.4; p.nc.intol=-1e-2; p=cont(p,35);% relax instab. tolerance  
p.nc.dsmax=0.025; p.nc.intol=0; p=cont(p,15); % tighten instab.tol, reduce max ds to capture fold 
%% 2nd pattern 
p=swibra(odir,'bpt2','c2',0.01);  p.nc.tol=1e-6; p.sw.bifcheck=0; 
p.nc.dsmax=0.3; p=cont(p,10); p.nc.dsmax=0.05; p=cont(p,20);
%% BD, max(u)
lws=2.5; lwu=1; ms=8; 
fnr=3; figure(3); clf; c=8; plotbra('h0','pt20',fnr,c,'lwst',lws,'lwun',lwu,'ms',ms,'bplab',[1 2]); 
plotbra('c1',fnr,c,'cl','b','lab',[10 30],'ms',0,'lwst',lws,'lwun',lwu,'lms',8); 
plotbra('c2','pt30',fnr,c,'cl','m','lab',10,'ms',0,'lwst',lws,'lwun',lwu); ylabel('max(u)'); 
%% BD , max(w)
fnr=3; figure(3); clf; c=10; plotbra('h0','pt20',fnr,c,'lwst',lws,'lwun',lwu,'ms',ms); 
plotbra('c1',fnr,c,'cl','b','lab',[10 30],'ms',0,'lwst',lws,'lwun',lwu); 
plotbra('c2','pt30',fnr,c,'cl','m','lab',10,'ms',0,'lwst',lws,'lwun',lwu); ylabel('max(w)'); 
%% BD, delta
fnr=3; figure(3); clf; c=7; plotbra('h0','pt20',fnr,c); 
plotbra('c1',fnr,c,'cl','b'); plotbra('c2','pt30',fnr,c,'cl','m'); ylabel('\lambda_m'); 
%% soln plots 
plotsol('c1','pt10'); figure(21); view(15,40); pause; plotsol('c1','pt30'); pause; plotsol('c2','pt10'); 
%%
p=loadp('c2','pt10'); p.plot.cview=[15,40]; p.plot.cut=[-2 0.001 -2]; 
plotsol(p); figure(21); %colorbar('Ticks',[1.6665,1.667]); 
colorbar('Ticks',[1.6362,1.6368]); 
colorbar('Ticks',[1.6462,1.6466]); 
