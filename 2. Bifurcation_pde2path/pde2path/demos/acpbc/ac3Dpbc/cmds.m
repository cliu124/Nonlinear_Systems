%% ac3Dpbc, but no transl.invariance,  demo to check basics
close all; keep pphome; rmdir('tr','s'); 
%% init 
p=[]; par=[1 -0.5 1 0]; lx=2*pi; ly=3*pi/2; lz=pi; per=[1 2 3]; 
p=acinit(p,lx,ly,lz,40,par,per); p.nc.nsteps=200; 
p.nc.ilam=2; p.nc.lammax=5; p.sol.ds=0.1; p.nc.dsmax=0.2; p=setfn(p,'tr'); 
%% use findbif to locate the BPs 
p=findbif(p,4);
%% switch to some bifurcating branches and continue
ind=[2 4]; 
for i=ind; 
    bp=['bpt' mat2str(i)]; out=['b' mat2str(i)]; 
    p=swibra('tr',bp,out,0.1); p=cont(p,10); 
end
%%
f=3; c=0; figure(f); clf;
plotbra('tr',f,c,'cl','c'); plotbra('b2','lab',10); 
plotbra('b4','pt10',f,c,'cl','r', 'lab',10); axis([-0.5 1 0 25]); 
%% solution plots 
v=[-13 70]; 
plotsol('b2','pt10',1,1,2); view(v); pause; % isosurf
plotsol('b2','pt10',1,1,1); view(v); pause; % slice 
plotsol('b4','pt10',1,1,2); view(v); % isosurf
%% fold-continuation, works as usual 
p=spcontini('b2','fpt1',3,'b2f');   % init fold continuation with par 3 new primary par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.01;                   % new stepsize in new primary parameter
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
%%
tic; p=cont(p,10); toc