%% cGL with NBC and time-periodic forcing, fixed period T, free nu 
% see nodalf.m and fofu.m for the forcing function, see cmds1b.m for free T 
close all; format compact; keep pphome; % clean up 
%% init, and continuation of trivial branch 
p=[]; lx=pi; nx=20;  
par=[-0.1; 1; 0.1; -1; 1; 0; 1;  0.5;    0.5];  
%    r,  nu,  mu, c3, c5, s,del, fo-ampl, forcing-cut 
p=cGLinit(p,lx,nx,par); dir='0a'; p=setfn(p,dir); % initialize  
p.usrlam=[]; p.nc.dsmax=0.2; p=cont(p,20); % run 
%% bifs at HBP1 and 2
figure(2); clf; aux=[]; aux.dlam=0; aux.tl=40; ds=0.1; dir='0a'; aux.freeT=0; 
for i=1:2
  hp=['hpt' mat2str(i)]; ndir=['sw' mat2str(i)];  
  p=hoswibra(dir,hp,ds,4,ndir,aux); p.nc.dsmax=0.4; 
  p=belon(p,2);  % Hopf border width is 3: pc, arclength, qh 
  p.hopf.ilam=2; p=cont(p,20); % free nu and go 
end 
%% use poiniguess
p=loadp('0a','pt0','sw3'); % load some steady point (for discr.data) 
nu=p.nu; p.u(nu+1)=0.5; % reset some parameter as desired
t=linspace(0,2*pi,40); ia1=0.0; ia2=0.1; % create guesses for IC 
%uv=ones(nu,1); % x-homo iguess, or choose the next line
x=getpte(p); x=x'; uv=[cos(x);cos(x)];  p.u(p.nu+1)=1; % an x-dependent guess
tl=length(t); u=zeros(nu,tl); 
for i=1:tl; u(:,i)=ia1+ia2*cos(t(i)+1)*uv; end
tausw=1; % 0=no tangent to iniguess (use u itself), 1: give tangent for u-vars, if 
%  desired adapt remaining vars afterwards. 
switch tausw; 
    case 0; aux=[]; aux.ds=0.1; p=poiniguess(p,t,u,aux); % 
    case 1; aux=[]; aux.ds=0.1; aux.tau=u; 
        p=poiniguess(p,t,u,aux); % variant with explicit tangent 
        p.hopf.tau(end-1)=-0.1; % reset some vars of tangent, here prim-cont param r; 
        % here dr<0 and ds>0 means that we aim to the left
end
p.hopf.freeT=0; p.hopf.ilam=2; % fix T and free nu 
p.sol.restart=0; % no restart (hence para=4, i.e., hopf.tau used) 
% if desired/needed, rather use an explicit tau in hoiniguess  to control first step: 
p.usrlam=[]; p.sw.verb=0; p=belon(p,2); p=cont(p,21); % go 
%% plot BD, max 
cmp=11; wnr=3; figure(wnr); clf; plotbra('sw1','pt20',wnr,cmp, 'lab', 20); 
plotbra('sw2','pt20',wnr,cmp, 'lab', 20,'cl','r'); 
plotbra('sw3','pt20',wnr,cmp, 'lab', 20,'cl','b'); 
xlabel('r'); ylabel('max(u)'); 
%% plot BD, nu
cmp=2; wnr=3; figure(wnr); clf; plotbra('sw1','pt20',wnr,cmp, 'lab', 20); 
plotbra('sw2','pt20',wnr,cmp, 'lab', 20,'cl','r'); 
plotbra('sw3','pt20',wnr,cmp, 'lab', 20,'cl','b'); 
xlabel('r'); ylabel('\nu'); 
%% soln plots
v=[10, 50]; 
hoplotf('sw1','pt20',1,1); figure(1); view(v); title('u_1 at n1/20'); pause 
hoplotf('sw2','pt20',1,1); figure(1); view(v); title('u_1 at n2/20'); pause 
hoplotf('sw3','pt20',1,1); figure(1); view(v); title('u_1 at n3/20'); 