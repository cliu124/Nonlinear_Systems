%% 10-node-tetras, long and slender bar, with snaking branch of BCC-front
keep pphome; close all;
%% init 
p=[];lx=pi/sqrt(2);ly=lx;lz=8*lx;nx=4; ny=round(ly/lx*nx); nz=round(lz/lx*nx); 
lam=-0.01;nu=1.5; par=[lam nu]; hofem.Kfn='K2'; hofem.sw=1; hofem.t2sinterpol=0; 
p=shinit(p,lx,ly,lz,nx,ny,nz,par,hofem); 
%% 
p=setfn(p,'0l'); p.sw.bifcheck=0; p.nc.neig=4; huclean(p); p.plot.pstyle=2; 
p.nc.eigref=-1; p.sw.verb=2; p.nc.bisecmax=3; p.sol.ds=0.1; p=cont(p,1);
%% 1st BP, double, spots and stripes 
aux.besw=1; aux.m=4;  p0=qswibra('0l','pt1',aux);
p0.nc.dsmax=0.025; p0=setbelilup(p0,1,1e-5,20,1e-3,100);
p0.nc.neig=4; p0.nc.eigref=-1; p0.sw.verb=2; p0.nc.dsmax=0.05;
p0.sw.bifcheck=0; p0.sw.foldcheck=0; p0.nc.bisecmax=4; 
%% hot balls 
huclean(p); p=seltau(p0,2,'b1l',2); p.sol.ds=-0.005; p.nc.dsmax=0.02; 
p.sw.bifcheck=1; tic; p=cont(p,15); toc
p.sw.bifcheck=0; p.nc.dsmax=0.05; tic; p=cont(p,10); toc 
%% 
p=loadp('b1l','bpt1'); p.plot.pstyle=3; p.sol.ds=-0.05; p.nc.dsmax=0.1; 
p.sw.bifcheck=0; p.sw.foldcheck=1; p=cont(p,50);
%% use gentau for tubes; 
p=gentau(p0,[0 0 1],'BCCtl'); pause; p=cont(p,20); 
%% swibra to localized BCCs 
p=swibra('b1l','bpt1','bf'); p.plot.pstyle=3; p.nc.dsmax=0.05; p=cont(p,50); 
%% plot BD 
f=3; c=4; figure(f); clf; plotbra('b1l','pt15',f,c,'lab',40,'cl','b'); 
plotbra('bf','pt160',f,c,'cl','r','lab',[40 140]); 
xlabel('\lambda'); ylabel('||u||_2'); 
axis([-0.35 0.02 0 0.45]); 
%% soln plot 
plotsol('bf','pt30',1,1,3); colormap hot; pause 
mypsol('bf','pt40'); pause; mypsol('bf','pt140'); 