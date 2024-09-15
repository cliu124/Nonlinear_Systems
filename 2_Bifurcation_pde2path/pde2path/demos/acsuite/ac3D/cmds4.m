%% AC3D on finer mesh to test different linear system solvers (LSSs) and Eval solvers 
close all; keep pphome; 
%% init, here exactly at 1st BP
p=[]; par=[1 0.4236 1 0]; lx=2*pi; ly=3*pi/2; lz=pi; 
nx=30; dir='tr1'; % nx=40; dir='tr2'; nx=50; dir='tr3'; % uncomment for the desired mesh
p=acinit(p,lx,ly,lz,nx,par); p.np, plotsol(p,1,1,1); p.fuha.spjac=@spjac; 
p.nc.ilam=2; p.nc.lammax=2; p.sw.spcalc=0; p.sol.ds=0.01; p.nc.dsmax=0.2; p=setfn(p,dir); p0=p; 
%% just do one (dummy step) to save pt 
p=p0; pp.sol.ds=1e-4; p=cont(p,1); % 1 step, only to generate tangent 
%% swibra to 1st branch, switch off: spectral comp, plotting, fold-detec, saving
p=swibra(dir,'pt1','dummy',0.05); p.sw.spcalc=0; p.plot.pmod=10; p.sw.foldcheck=0; 
p.file.smod=0; p0=p; 
%% Now test LSSs, first standard \: nx=30: 19s, nx=40: 85s, nx=50: 225
p=p0; tic; p=cont(p,30); toc 
%% bel with \: nx=30: 18s,  nx=40: 80s, nx=50: 210s 
p=p0; bw=0; beltol=1e-3; belmaxit=5; p=setbel(p,bw,beltol,belmaxit,@lss); 
tic; p=cont(p,30); toc 
%% lsslu, lsstol=1e-3; speed essentially like \;  lsstol=1e-2  ->  too small ds; 
p=p0; p.fuha.lss=@lsslu; tic; p=cont(p,30); toc 
%% bel with AMG: nx=30: 11s, nx=40: 30s, ilun=2: 38s, nx=50: 85s
p=p0; bw=0; beltol=1e-3; belmaxit=5; droptol=1e-3; amgmaxit=50;  % param. for blssbel,  
p=setbelilup(p,bw,beltol,belmaxit,droptol,amgmaxit); % lssAMG as inner LSS
p.sw.verb=3; p.ilup.ilun=0;  % report timing of lssAMG, no forced prec updates 
%p.ilup.ilun=2; p=setfn(p,'p1c2'); % uncomment this line for prec update each ilun-th-steps
tic; p=cont(p,30); toc 
%% \, with spcalc by standard eigs. nx=30: 31s, nx=40: 123s, nx=50: 350s
p0.sw.spcalc=1; p=p0;  tic; p=cont(p,30); toc 
%% bel with AMG + eigs: nx=30: 23s, nx=40: 80s, nx=50: 220s 
p=p0; p=setbelilup(p,bw,beltol,belmaxit,droptol,amgmaxit); % lssAMG as inner LSS
p.sw.verb=3;  p.sw.eigssol=0; tic; p=cont(p,30); toc 
%% bel with AMG + ilueigs, default: nx=30: 80s, nx=40: 213s, nx=50: 440s 
% (cause eigs becomes slow when prec is poor)
% with ilun=2: nx=30: 35s; nx=40: 50s,  nx=50: 160s
p=p0; p=setbelilup(p,bw,beltol,belmaxit,droptol,amgmaxit); % lssAMG as inner LSS
p.sw.verb=2; p.sw.spcalc=1; p.sw.eigssol=2; p.sw.bifcheck=2; 
p.ilup.ilun=2;  % uncomment this line for prec-update each ilun-th-steps
tic; p=cont(p,30); toc 
%% check LSS for fold-continuation; 
% swibra to 1st branch, switch off spectral comp, plotting, fold-detec, saving
dir='tr1'; odir='fc1'; %dir='tr2'; odir='fc2'; dir='tr3'; odir='fc3'; 
p=swibra(dir,'pt1',odir,0.05); p.sw.spcalc=1; p.plot.pmod=10; p.sw.foldcheck=1; 
p.file.smod=0; p0=p; 
%% AMG: nx=30: 11s, nx=40: 33s, ilun=2: 38s, nx=50: 86s
p=p0; bw=0; beltol=1e-3; belmaxit=5; droptol=1e-3; amgmaxit=50;  % param. for blssbel,  
p=setbelilup(p,bw,beltol,belmaxit,droptol,amgmaxit); % lssAMG as inner LSS
p.sw.verb=1; p.sw.ilumod=0; p.sw.foldcheck=1; 
%p.ilup.ilun=2; p=setfn(p,'p1c2'); % uncomment this line for ilu-update each ilun-th-steps
tic; p=cont(p,30); toc 
%% fold-continuation. 
p=spcontini(odir,'fpt1',3,'dummy');   % init fold continuation with par 3 new primary par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new parameter for plotting
p.sol.ds=-0.1;                   % new stepsize in new primary parameter
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
p.sw.spcalc=0; p.nc.lammin=0.5; p.sw.bifcheck=0; p.sw.foldcheck=0; 
p0=p; 
%% bel with AMG:  nx=30: 5s, nx=40: 12s 
p=p0; p.sw.verb=3; tic; p=cont(p,15); toc
%% back to \: nx=30: 243s, nx=40: 1700s, bel with \ no better
p=p0; p.fuha.lss=@lss; p.fuha.blss=@lss; tic; p=cont(p,10); toc
