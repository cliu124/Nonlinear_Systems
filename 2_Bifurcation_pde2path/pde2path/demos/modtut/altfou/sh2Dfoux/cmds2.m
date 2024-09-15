%% SH 2D via dct, larger scale (slow, rather see sh2dfou and sh2dmfree !) 
close all; keep pphome; 
%% init 
p=[]; par=[-0.02 2 -1]; % lam,quad,cubic, 
lx=12*pi; ly=sqrt(3)*lx/3; nx=100; ny=round(nx*ly/lx); 
p=shinit(p,lx,ly,nx,ny,par,Lfn); 
p=setfn(p,'tr2'); p.sw.verb=2; p.ps=2; p.sw.bifcheck=0; % switch off bifcheck for speed  
p.sw.spcalc=1; p.nc.eigref=-0.1; 
tic; p=cont(p,2); toc
%% 1BP double 
aux=[]; aux.soltol=1e-10; aux.m=2;  aux.isotol=1e-12; 
p0=qswibra('tr2','pt2',aux); p0.nc.tol=1e-6;% cswibra, then reset some parameters
%% select tangent and cont. 
p=seltau(p0,2,'hexl',2); p.sol.ds=-0.01; p.nc.dsmax=0.05; p.nc.mu2=0.005; 
p.sw.bifcheck=2; p.nc.bisecmax=3; p.sw.eigssol=0; p=cont(p,5); pause
p.sw.bifcheck=0; p.nc.tol=1e-8; p=cont(p,20); 
%% stripes via gentau 
p=gentau(p0,[0 1]);  p=setfn(p,'b1l'); p.nc.tol=1e-6; p.sol.ds=-0.1;  p.nc.dsmax=0.2; 
p.nc.dsmin=0.01; p.sw.bifcheck=0; p=cont(p,1); p.nc.tol=1e-8; p=cont(p,20); 
%% hexagon-front 
p=swibra('b2l','bpt1','hf'); p.sw.bifcheck=0; 
p.nc.dsmax=0.2; p.nc.dsmin=0.001; p.sw.foldcheck=0; 
%%
ilu=0; % ilu works but not very efficient here 
if ilu; p.sw.verb=3; p=setbelilup(p,0,1e-4,5,1e-6,200); 
else p.fuha.innerlss=@lss; end %p=setbel(p,0,1e-6,5,@lsscg); end 
p.nc.neig=1; p.nc.eigref=-2; % just check stab. 
tic; p=cont(p,10); toc 
%% BD plot 
f=3; c=5; figure(f); clf;
plotbra('tr2',f,c,'cl','k','lsw',0); 
plotbra('hexl',f,c,'cl','b','lab',20); 
plotbra('hf',f,c,'cl','r','lab',[40 80],'lp',128); 
ylabel('||u||_2'); 
%% soln plot 
p=loadp('hexl','pt45'); p.ps=3; plotsol(p);
%%
plotsol('hf','pt40'); pause; plotsol('hf','pt80'); 