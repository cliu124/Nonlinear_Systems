%% script for Hopf bif in pollution model Wirl2000, here ODE case 
close all; keep pphome 
%% C1: init and continue trivial branch 
p=[]; par=[0.5 1 0.2 0 300]; % [del, pr, beta, a, ga]; 
p=pollinit(p,par); p=setfn(p,'FSS'); screenlayout(p); p.file.smod=2; 
p.sw.bifcheck=2; p.sw.verb=2; p.nc.mu2=1e-3; % accuracy of Hopf detection 
p.nc.ilam=1; p.sol.ds=0.01; p.nc.dsmax=0.01; p.usrlam=0.54:0.01:0.6; p0=p; 
%%
p=cont(p,20); % cont of FSS
%% C2: cont of Hopf branches 
para=4; ds=0.5; dsmax=1; xi=1e-2; figure(2); clf; aux=[]; aux.tl=80;  
p=hoswibra('FSS','hpt1',ds,para,'h1',aux); nsteps=20; pause
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.hopf.y0dsw=0; 
p.file.smod=1; p.hopf.flcheck=0; p.sw.verb=1; 
tic; p=cont(p,nsteps);  toc
%% C3: plot the BD, data on branch is 
% par(1..5), T, min, max, |u(.,.)|_L^2, J_c(CSS) resp J_c(u_H), J_c(u_H(.+T/2))
%            6                            10                         11         
figure(3); clf; pcmp=10; % use var for plot-cmp for quick changing 
plotbra('FSS','pt20',3,pcmp,'cl','k','lsw',1,'ms',0); 
plotbra('h1','pt20',3,pcmp,'cl','r','lab',[13]); 
plotbra('h1','pt20',3,11,'cl','r','tyun','--','lsw',0); 
axis([0.5 0.6 -0.05 0.35]); xlabel('\rho'); ylabel('J'); 
%% C7: floqps works reasonably 
[muv1, muv2,ind,h]=floqpsap('h1','pt13');
%%
fprintf('%g %g %g\n',muv1(end-1),muv1(end),muv2(1)); 
%%
ps=loadp('FSS','pt13'); plotsol(ps,1,1,1); 