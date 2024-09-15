%% script for Hopf bif in pollution model Wirl2000, here with diffusion 
close all; keep pphome 
%% C1: init and continue trivial branch 
p=[]; lx=pi/2; nx=20; par=[0.5 1 0.2 0 300]; % [del, pr, beta, a, ga]; 
p=pollinit(p,lx,nx,par); p=setfn(p,'FSS'); screenlayout(p); p.file.smod=2; 
p=initwn(p,2,1); p=initeig(p,1); p.nc.neig=[10 10]; % find guess for omega_1 
p.sw.bifcheck=2; p.sw.verb=2; p.nc.mu2=1e-3; % accuracy of Hopf detection 
p.nc.ilam=1; p.sol.ds=0.01; p.nc.dsmax=0.01; p.usrlam=0.54:0.01:0.6; p0=p; 
%%
p=p0; p=cont(p,10); % cont of FSS
%% C2: cont of Hopf branches 
para=4; ds=0.5; dsmax=1; xi=1e-2; figure(2); clf; aux=[]; aux.tl=40;  
for j=1:2
 switch j 
     case 1; p=hoswibra('FSS','hpt1',ds,para,'h1',aux); nsteps=20; pause
     case 2; p=hoswibra('FSS','hpt3',ds,para,'h2',aux); nsteps=20; pause
 end
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.hopf.y0dsw=0; 
p.file.smod=1; p.hopf.flcheck=0; % use floqps for multipliers 
%p.usrlam=0.54:0.01:0.56; 
p.sw.verb=1; 
tic; p=cont(p,nsteps);  toc
end 
%% C3: plot the BD, data on branch is 
% par(1..5), T, min, max, |u(.,.)|_L^2, J_c(CSS) resp J_c(u_H), J_c(u_H(.+T/2))
%            6                            10                         11         
figure(3); clf; pcmp=10; % use var for plot-cmp for quick changing 
plotbra('FSS','pt20',3,pcmp,'cl','k','lsw',0,'ms',0); 
plotbra('h2','pt25',3,pcmp,'cl','r','lab',[17]); 
plotbra('h2','pt25',3,11,'cl','r','tyun','--','lsw',0); 
plotbra('h1','pt20',3,pcmp,'cl','b', 'lab',8); 
axis([0.5 0.6 -0.05 0.35]); xlabel('\rho'); ylabel('J'); 
%% C4: solution plots for h1
hoplotf('h1','pt8',1,1); figure(1); title(''); %title('h1/pt8'); 
xlabel('x'); ylabel('t');  zlabel('emissions'); pause
hoplotf('h1','pt8',1,2); figure(1); title(''); %title('h1/pt8'); 
xlabel('x'); ylabel('t');  zlabel('stock'); pause
hoplotf('h1','pt8',1,5); figure(1); title(''); %title('h1/pt8'); 
xlabel('x'); ylabel('t');  zlabel('J_c'); pause
hoplotf('h1','pt8',1,6); figure(1); title(''); %title('h1/pt8'); 
xlabel('x'); ylabel('t');  zlabel('control'); 
%% C5: solution plots for h2
hoplotf('h2','pt17',1,1); figure(1); title('h2/pt17'); 
xlabel('x'); ylabel('t');  zlabel('emissions'); pause
hoplotf('h2','pt17',1,2); figure(1); title('h2/pt17'); 
xlabel('x'); ylabel('t');  zlabel('stock'); pause
hoplotf('h2','pt17',1,5); figure(1); title('h2/pt17'); 
xlabel('x'); ylabel('t');  zlabel('J_c'); pause
hoplotf('h2','pt17',1,6); figure(1); title('h2/pt17'); 
xlabel('x'); ylabel('t');  zlabel('control'); 
%% C6: illustrate full spectrum! 
p=loadp('FSS','pt2');
p.nc.neig=p.nu; r=resi(p,p.u); Gu=getGupde(p,p.u,r); 
[ineg,muv]=spcalc(Gu,p); 
%% C7: floqps works reasonably 
[muv1, muv2,ind,h]=floqpsap('h2','pt17'); fprintf('d(u_H)=%i\n',ind);
fprintf('%g %g %g\n',muv1(end-1),muv1(end),muv2(1)); 