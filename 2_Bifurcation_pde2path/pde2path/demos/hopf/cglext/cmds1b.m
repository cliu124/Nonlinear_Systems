%% cGL with NBC and time-periodic forcing, like cmds1, but with free T 
% and fixed nu 
close all; format compact; keep pphome; % clean up 
%% bifs at HBP1 & 2
figure(2); clf; aux=[]; aux.dlam=0; aux.tl=40; nsteps=20; ds=0.1; dir='0a'; 
for i=1:2
  hp=['hpt' mat2str(i)]; ndir=['s' mat2str(i)]; 
  p=hoswibra(dir,hp,ds,4,ndir,aux); p.nc.dsmax=0.4;
  p=belon(p,2); p.sw.verb=0; p=cont(p,nsteps); % go 
end 
%% use poiniguess
p=loadp('0a','pt0','s3'); % load some steady point (for discr.data) 
huclean(p); nu=p.nu; p.u(nu+1)=0.5; % reset some parameter as desired
T=2*pi; t=linspace(0,T,40); % create guesses for IC (incl. period T) 
ia1=0.0; ia2=0.1; % amplitudes for Iguess; ia2>0 needed here!!!
x=getpte(p); x=x'; uv=[cos(x);cos(x)];  p.u(p.nu+1)=1; % an x-dependent guess
tl=length(t); u=zeros(nu,tl); 
for i=1:tl; u(:,i)=ia1+ia2*cos(t(i)+1)*uv; end
aux=[]; aux.ds=0.1; p=poiniguess(p,t,u,aux); p.sol.restart=0; 
p.usrlam=[]; p.sw.verb=0; p.sw.bifcheck=0; p=cont(p,21); % go 
%% plot BD, max 
cmp=11; wnr=3; figure(wnr); clf; plotbra('s1','pt20',wnr,cmp, 'lab', 20); 
plotbra('s2','pt20',wnr,cmp, 'lab', 20,'cl','r'); 
plotbra('s3','pt20',wnr,cmp, 'lab', 20,'cl','b'); 
xlabel('r'); ylabel('max(u)'); 
%% plot BD, T
cmp=10; wnr=3; figure(wnr); clf; plotbra('s1','pt20',wnr,cmp, 'lab', 20,'fp',1); 
plotbra('s2','pt20',wnr,cmp, 'lab', 20,'cl','r','fp',1); 
plotbra('s3','pt20',wnr,cmp, 'lab', 20,'cl','b'); 
xlabel('r'); ylabel('T'); 
%% soln plots
v=[10, 50]; 
hoplotf('s1','pt20',1,1); figure(1); view(v); title('u_1 at s1/20'); pause 
hoplotf('s2','pt20',1,1); figure(1); view(v); title('u_1 at s2/20'); pause 
hoplotf('s3','pt20',1,1); figure(1); view(v); title('u_1 at n3/20'); 