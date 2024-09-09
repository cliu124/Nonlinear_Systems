% cgl, 2D, DBC 
close all; keep pphome; 
%% init
ndim=2; dir='hom2d'; p=[]; lx=pi; nx=40; par=[0.75; 1; 0.1; -1; 1]; 
p=cGLinit(p,lx,nx,par,ndim); p=setfn(p,dir); p.nc.mu1=1; 
p.sw.bifcheck=2; p=cont(p,20);
%% first 2 branches 
para=4; ds=0.1; aux=[]; aux.tl=30; figure(2); nsteps=15; clf; 
for bnr=1:2
switch bnr
 case 1; p=hoswibra('hom2d','hpt1',ds,para,'2db1',aux); 
 case 2; p=hoswibra('hom2d','hpt2',ds,para,'2db2',aux); 
end 
p.hopf.jac=1; p.nc.dsmax=0.4; p.hopf.xi=1e-3; p.sw.verb=2; p.hopf.flcheck=0; 
AMG=0;  % set AMG=1 if ilupack available
if ~AMG; p=setbel(p,1,1e-4,20,@lss); % use BEL without ilupack 
else  % use AMG with or without bel; AMG seems indifferent to borders! 
    p=setilup(p,1e-3,100); p.fuha.lss=@lss; p.fuha.blss=@lssAMG;
end 
tic; p=cont(p,nsteps);  toc
end 
%% BD, L^2
bpcmp=9; figure(3); clf; plotbra('hom2d',3,bpcmp,'cl','k','lsw',0); 
plotbra('2db1',3,bpcmp,'cl','r'); 
plotbra('2db2',3,bpcmp,'lab',8,'cl','b'); 
axis([1 2.6 0 1]); xlabel('r'); ylabel('||u||_*'); 
%% BD, T 
bpcmp=6; figure(3); clf; plotbra('2db1',3,bpcmp,'cl','k'); 
plotbra('2db2',3,bpcmp,'lab',8,'cl','b'); 
axis([1 2.3  6.33 7.5]); xlabel('r'); ylabel('T');
%% soln plot
aux=[];  aux.pstyle=3; aux.xtics=[-2 2]; aux.ytics=[-1 1]; 
hoplotf('2db2','pt8',1,1,aux); 
%% stab-cmds 
p=loadp('2db2','pt8'); hoplot(p,1,1); dir='stab2d'; 
p.u(1:p.nu)=p.hopf.y(1:p.nu,1); u0=p.u(1:p.nu); p=setfn(p,dir);
ts=[]; t0=0; npp=40; pmod=50; smod=5; tsmod=1; nc=0; 
%% repeat if necessary ...  
nt=400;
[p,t0,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,@nodalf,1); 
figure(4); clf; plot(ts(1,:), ts(2,:)); % ts-plot, point-vals
figure(5); clf; plot(ts(1,:), ts(3,:)); % difference in norm
set(gca,'FontSize',p.plot.fs); axis tight; axis([0 10*p.hopf.T 0 3]); 
xlabel('t'); ylabel('||u(t)-u_0||_{\infty}'); 
%% soln plot
dir='stab2d'; si=0; incr=80; nt=1*npp/smod; wnr=2; cmp=1; pstyle=2; 
ind=si:incr:5*incr; lay=[2 3]; notic=1; 
tintplot2d(dir,ind,wnr,cmp,pstyle,lay,notic);
%% movie
p=loadp('2db2','pt8'); aux.xtics=[-2 2]; aux.ytics=[-1 1]; 
mov=homov2d(p,1,1,aux); mymov2avi(mov,'mcGL2D');