%% script for Hopf bif in pollution model Wirl2000, with diffusion, outdated 
% cause extended to include control towards POs in ocdemos/pollution
close all; keep pphome 
%% C1: init and continue trivial branch 
p=[]; lx=pi/2; nx=40; par=[0.5 1 0.2 0 300]; % [del, pr, beta, a, ga]; 
p=pollinit(p,lx,nx,par); p=setfn(p,'FSS'); screenlayout(p); p.file.smod=2; 
p=initwn(p,2,1); p=initeig(p); p.nc.neig=[5 5]; % find guess for omega_1 
p.sw.bifcheck=2; p.sw.verb=2; p.nc.mu2=1e-3; % accuracy of Hopf detection 
p.nc.ilam=1; p.sol.ds=0.01; p.nc.dsmax=0.01; p=cont(p,20); % cont of FSS
%% C2: cont of Hopf branches 
para=4; ds=0.5; dsmax=1; xi=1e-2; figure(2); clf; aux=[]; aux.tl=25;  
for j=1:2
 switch j 
     case 1; p=hoswibra('FSS','hpt1',ds,para,'h1',aux); nsteps=15; 
     case 2; p=hoswibra('FSS','hpt2',ds,para,'h2',aux); nsteps=25; 
 end
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.hopf.y0dsw=0; 
p.sw.verb=1; p.file.smod=1; p.hopf.flcheck=2; % use floqps for multipliers 
p.usrlam=[0.5 0.6 0.7]; tic; p=cont(p,nsteps);  toc
end 
%% C3: plot the BD, data on branch is 
% par(1..5), T, min, max, |u(.,.)|_L^2, J_c(CSS) resp J_c(u_H), J_c(u_H(.+T/2))
%            6                            10                         11         
figure(3); clf; pcmp=10; % use var for plot-cmp for quick changing 
plotbra('FSS','pt20',3,pcmp,'cl','k'); 
plotbra('h2','pt25',3,pcmp,'cl','r','lab',[4,14]); 
plotbra('h2','pt25',3,11,'cl','r','tyun','--','lsw',0); 
plotbra('h1','pt16',3,pcmp,'cl','b', 'lab',4); 
axis([0.5 0.61 -0.05 0.35]); xlabel('\rho'); ylabel('J'); 
%% C4: solution plots for h1
hoplotf('h1','pt4',1,1); figure(1); title('h1/pt4'); 
xlabel('x'); ylabel('t');  zlabel('emissions'); pause
hoplotf('h1','pt4',1,2); figure(1); title('h1/pt4'); 
xlabel('x'); ylabel('t');  zlabel('stock'); pause
hoplotf('h1','pt4',1,5); figure(1); title('h1/pt4'); 
xlabel('x'); ylabel('t');  zlabel('J_c'); pause
hoplotf('h1','pt4',1,6); figure(1); title('h1/pt4'); 
xlabel('x'); ylabel('t');  zlabel('control'); 
%% C5: solution plots for h2
hoplotf('h2','pt14',1,1); figure(1); title('h2/pt14'); 
xlabel('x'); ylabel('t');  zlabel('emissions'); pause
hoplotf('h2','pt14',1,2); figure(1); title('h2/pt14'); 
xlabel('x'); ylabel('t');  zlabel('stock'); pause
hoplotf('h2','pt14',1,5); figure(1); title('h2/pt14'); 
xlabel('x'); ylabel('t');  zlabel('J_c'); pause
hoplotf('h2','pt14',1,6); figure(1); title('h2/pt14'); 
xlabel('x'); ylabel('t');  zlabel('control'); 
%% C6: illustrate full spectrum! 
p=loadp('FSS','pt0');
p.nc.neig=p.nu; r=resi(p,p.u); Gu=getGupde(p,p.u,r); 
[ineg,muv]=spcalc(Gu,p); 
%% C7: floqps works reasonably 
[muv1, muv2,ind,h]=floqpsap('h2','pt4'); fprintf('d(u_H)=%i\n',ind); 