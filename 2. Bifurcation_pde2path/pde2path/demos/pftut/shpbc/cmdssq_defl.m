%% SH on square via deflation 
close all; keep pphome; 
%% init at lam=1  (rather strongly supercritical)
p=[]; lx=2*pi; nx=round(6*lx); ly=lx; dir='sq2'; ndim=2; lam=1; nu=0; 
par=[lam; nu; 0; 0];  % s1,s2 for PC 
sw.ref=0; sw.sym=2; p=shinit(p,nx,lx,ly,ndim,par,sw); huclean(p); p=setfn(p,dir); 
p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.02; p.np
p.file.smod=10; p.nc.tol=1e-8; p.nc.lammax=2; 
%% prepare deflation 
x=getpte(p); x=x'; x=[x; x]; x=p.mat.drop*x; p.u(1:p.nu)=zeros(p.nu,1); 
x2=max(x(:,1)); y2=max(x(:,2)); u0=p.u; 
%p.defl.al=0.1; p.defl.al2=1e3; p.defl.nsw=2; 
al3=0; % al3=0: use SMW-solver,  al3>0: use lumping (with factor al3) 
global p2pglob; p=deflinit(p,al3); 
%% run deflated Newton with iguesses from random F-coeffis 
nu=p.nu; amp=1; p.nc.imax=20; p.nc.tol=1e-4; 
p.defl.nd=1; % number of known solns 
p.defl.nsw=2; ng=10; % norm-switch (||u||_L^2) and # of guesses 
p.sw.newt=0; % newton-switch: 0/1=nloop, 2=NLEQ1, 3=fsolve
ww=[[0.25;0] [0; 0.25] [0.25; 0.25]  [0.5;0] [0; 0.5] [0.5; 0.5] ... % 1-6
    [0.5;0.25] [0.25; 0.5] [1;0] [0;1] [1;1] [1; 0.25] [0.25;1]  ... % 7-13    
    [1; 0.5] [0.5;1] [1.25;0] [0;1.25] [1.25;0.25] [1.25;1.25]];     % 14-19
wnr=sum(ww.^2,1); 
wi=10:18; %indices of wave vectors, other choices: % wi=1:2:3; % long, wi=12:14; % special 
for ni=1:ng % deflation loop 
  fc=rand(length(ww),1)-0.5; ug=zeros(nu,1); 
  for i=wi;  % get ICs from random Fourier
    sx=x2*(rand(1)-0.5); sy=y2*(rand(1)-0.5); 
    for j=1:nu/2       
        ug(j)=ug(j)+fc(i)*cos((x(j,:)-[sx, sy])*ww(:,i));
        ug(j+nu/2)=ug(j+nu/2)-wnr(i)*fc(i)*cos((x(j,:)-[sx, sy])*ww(:,i));
    end
  end 
  u0(1:p.nu)=amp*ug; p.u=u0; plotsol(p); %pause
  fprintf('i=%i ',ni); ok=1; 
  while ok>0; ndo=p.defl.nd; 
    [u, p]=deflsol(p,u0);% p.u=u; plotsol(p); 
    if p.defl.nd==ndo; ok=0; end
  end
end
%%
p.defl.nd, for i=2:p.defl.nd; p.u=p.defl.u(:,i); plotsol(p); nolabti; title(mat2str(i)); pause; end
%% save defl soln in p1 and run cont on selected ones 
p.nc.dsmax=0.1; p1=p; % save p in p1 to start cont. runs
p1.nc.lammax=2; p1.sw.spcalc=1; p.usrlam=[0 1 2]; 
p=postdefl(p1,2,'d2a',0.1); p=qxyon(p); p.sw.spcalc=1; p=cont(p,20); 
p=postdefl(p1,2,'d2b',-0.1); p=qxyon(p); p=cont(p,40);
%%
%p=postdefl(p1,3,'j3a',0.1); p=qxyon(p); p=cont(p,20); 
p=postdefl(p1,3,'j3b',-0.1); p=qxyon(p); p=cont(p,40); 
%% 
f=3; c=7; figure(f); clf;
plotbra('d2a',f,c,'cl','r','lsw',0); 
plotbra('d2b',f,c,'cl','r','lab',38);
%plotbra('j3a',f,c,'cl','b','lsw',0); 
%plotbra('j3b',f,c,'cl','b','lab',29);
%axis([-0.2 0.51 0.1 3.5]); 
xlabel('\lambda'); ylabel('min(u)'); 
%%
plotsol('d2a','pt30'); nolabti;

