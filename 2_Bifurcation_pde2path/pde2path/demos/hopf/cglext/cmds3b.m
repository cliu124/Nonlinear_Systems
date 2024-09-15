%% cGL with additive x-dep.forcing (fofu2), freeT=0, pc=0; 
close all; format compact; keep pphome; % clean up 
%% just init, i.e., generate x-mesh and data structures, then save  
p=[]; lx=pi; nx=20;  
par=[-0.1; 1;  0.5; 0.2;  0;  0; 1];  
%    r,   nu,  mu,  c3,  c5,   s, del=x-wave-nr of forcing 
p=cGLinit(p,lx,nx,par); dir='0c'; p=setfn(p,dir); % initialize  
p.fuha.sG=@sG2; p.fuha.sGjac=@sGjac2; p.usrlam=[]; p.nc.dsmax=0.2; 
stansavefu(p);  p.nc.ilam=1; p=cont(p,20); 
%% bifs at HBP1 and 2
figure(2); clf; aux=[]; aux.dlam=0; aux.tl=40; ds=0.1; dir='0c'; aux.freeT=0; 
for i=1:2
  hp=['hpt' mat2str(i)]; ndir=['csw' mat2str(i)];  
  p=hoswibra(dir,hp,ds,4,ndir,aux); p.nc.dsmax=0.4; 
  p=belon(p,2); p.hopf.flcheck=0;  % Hopf border width is 3: pc, arclength, qh 
  p.hopf.ilam=2; p=cont(p,10); % free nu and go 
end 
%% switch to cont in c5, no PC, fix nu 
p=hoswiparf('csw1','pt10','csw1c',5,0.01); p.hopf.pc=0;  p.hopf.ilam=[]; 
%p.hopf.tau=0*p.hopf.tau; p.hopf.tau(end-1:end)=-0.1; 
p.nc.dsmax=0.3; p=cont(p,21); 
%% switch to cont in c5, no PC, fix nu 
p=hoswiparf('csw2','pt10','csw2c',5,0.01); p.hopf.tau=0*p.hopf.tau; 
%%
p.hopf.tau(end-1:end)=0.1; p.hopf.pc=0; p.hopf.ilam=[]; 
p.nc.dsmax=0.3; p=cont(p,3); 
%% plot BD, max (9), nu (2) etc  
cmp=9; %cmp=2; 
wnr=3; figure(wnr); clf; lli1=[10 20]; lli2=[10 30]; 
plotbra('csw1c','pt20',wnr,cmp,'cl','b','lab',lli1,'fp',1); 
plotbra('csw2c','pt30',wnr,cmp,'cl','r','lab',lli2,'fp',1); 
xlabel('c_5'); ylabel('max(u)'); 
%% soln plots
v=[40, 40]; cmp=1; 
hoplotf('csw1c','pt10',1,cmp); figure(1); view(v); title('u_1 at c1/10'); pause 
hoplotf('csw1c','pt20',1,cmp); figure(1); view(v); title('u_1 at c1/20'); pause; 
hoplotf('csw2c','pt10',1,cmp); figure(1); view(v); title('u_1 at c2/10'); pause 
hoplotf('csw2c','pt30',1,cmp); figure(1); view(v); title('u_1 at c2/30'); 
