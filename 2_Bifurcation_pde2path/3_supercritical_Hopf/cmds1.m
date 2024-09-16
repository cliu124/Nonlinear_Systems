%% demo for supercritical Hopf bifurcation on page 73 of Khalil (2002)
% This is a demo for 
% \dot{x}_1=x1*[\mu-(x1^2+x2^2)]-x2
% \dot{x}_2=x2*[\mu-(x1^2+x2^2)]+x1

close all; keep pphome; 
%% cell 1: init and cont of trivial branch 
p=[]; par=[-2]; % here par(4)=coefficient of x-dependent terms 
p=init(p,par); p=setfn(p,'tr'); p=cont(p);
%% cell 2: switch to first 3 bifurcating branches and continue
para=4; ds=0.1; figure(2); clf; aux=[]; 
p=hoswibra('tr','hpt1',ds,para,'hpt1',aux); nsteps=20;
p.hopf.nfloq=2;
p.hopf.jac=1; p.nc.dsmax=0.5; p.hopf.xi=0.05; p.file.smod=5; p.sw.verb=2; 
p.hopf.flcheck=1; % switch for Floquet-comp: 0: off, 1:floq, 2: floqps 
p.sw.bifcheck=1;  % switch for bifurcation detection: 0:off, 1:on  
p=hocont(p,nsteps);
plotbra('tr'); plotbra('hpt1');