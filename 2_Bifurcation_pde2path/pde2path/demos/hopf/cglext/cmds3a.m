%% cGL with additive forcing (fofu2), freeT=0, pc=0; 
close all; format compact; keep pphome; % clean up 
%% just init, i.e., generate x-mesh and data structures, then save  
p=[]; lx=pi; nx=20;  
par=[0.5; 1;  0.5; 0.2;  0.5;  0; 0];  
%    r,   nu,  mu, c3,  c5,   s, del, 
p=cGLinit(p,lx,nx,par); dir='1a'; p=setfn(p,dir); % initialize  
p.fuha.sG=@sG2; p.fuha.sGjac=@sGjac2; p.usrlam=[]; p.nc.dsmax=0.2; 
stansavefu(p);  
%% poiniguess 
p=loadp('1a','pt0','a1'); huclean(p); nu=p.nu; 
T=2*pi; t=linspace(0,T,40); % create guesses for IC (incl. period T) 
ia1=0; ia2=0.5; uv=ones(nu,1); % x-homo iguess, 
tl=length(t); u=zeros(nu,tl); 
for i=1:tl; u(:,i)=ia1+ia2*sin(t(i)+2)*uv; end % any shift gives soln with same phase 
aux=[]; aux.ds=0.1; p=poiniguess(p,t,u,aux); % 
p.sol.restart=0; p.sw.bifcheck=1; p.hopf.flcheck=1; p0=p; 
% find initial soln, with p.hopf.pc=0 (non-autonomous, no PC, fixed T) 
p=p0; p.sw.verb=2; p.hopf.pc=0; p.hopf.freeT=0; 
[y,T,lam,res,iter,A,p]=honloop(p,p.hopf.y,p.hopf.T,p.hopf.lam); 
p.hopf.y=y; hoplot(p,11,1,1); p0=p; 
%% cont 
p=p0; p.nc.tol=1e-6; huclean(p);  p.usrlam=1; p.nc.dsmax=0.5; p=cont(p,10); 
%% switch to cont in c5 
p=hoswiparf('a1','pt8','a2a',5,0.1); p.nc.dsmax=0.2; p.sw.verb=0; p=cont(p,200); 
%% 2ndary POs, bifurcating from spat.homogen. one
p=poswibra('a2a','bpt2','s2',2); p=cont(p,20);
%% plot BD, max 
cmp=9; wnr=3; figure(wnr); clf; plotbra('a2a','pt200',wnr,cmp, 'lab', [0 70 180]); 
plotbra('s2','pt20',wnr,cmp,'cl','m', 'lab', 5); 
xlabel('c_5'); ylabel('max(u)'); 
%% soln plots
v=[40, 40]; 
hoplotf('a2a','pt0',1,1); figure(1); view(v); title('u_1 at a/0'); pause 
hoplotf('a2a','pt70',1,1); title('u_1(-\pi,t) at a/70'); pause 
hoplotf('a2a','pt180',1,1); title('u_1(-\pi,t) at a/180'); xlabel('t'); pause 
hoplotf('s2','pt5',1,1); figure(1); view(v); title('u_1 at s2/5'); 
%% change phase, falls back to org phase
p=loadp('a1','pt7'); 
y0=p.hopf.y; y1=circshift(y0(:,1:end-1),15,2); y1=[y1 y1(:,1)]; 
y0d=sety0dot(p,y,par,T); p.hopf.y0dsw=1; p.hopf.y0d=y0d; %p.hopf.pcfac=1e-3; 
p.hopf.y=y1; hoplot(p,10,1,1); 
[y,T,lam,res,iter,A,p]=honloop(p,p.hopf.y,p.hopf.T,p.hopf.lam); 
p.hopf.y=y; hoplot(p,11,2,1); 
