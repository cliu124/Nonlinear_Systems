%% demo to compare \ and lssbel with lssAMG as inner solver
close all; keep pphome; 
%% C1: init (generic), then specific settings (could also be set in init)
p=[]; par=[1 3 1 0]; lx=pi/2; ly=pi/2; lz=pi/2; nx=10; % choose different nx here
p=acinit(p,lx,ly,lz,nx,par); p.np, plotsol(p,1,1,1); p.fuha.spjac=@spjac; 
p.nc.ilam=2; p.nc.lammax=5; p.sol.ds=0.1; p.nc.dsmax=1; 
dir=mat2str(nx); p=setfn(p,dir); p.sw.verb=3; % outdir, and high verbosity 
p.sw.spcalc=0; p.sw.bifcheck=0; p.sw.foldcheck=0; % switch off spectral stuff 
bw=0; beltol=1e-3; belmaxit=5; droptol=1e-3; amgmaxit=50;  % param. for blssbel,  
p=setbelilup(p,bw,beltol,belmaxit,droptol,amgmaxit); % lssAMG as inner LSS
p=cont(p,1); % 1 step, only to generate tangent 
%% C2: localize BPs by setting lam to near a (here known) BP and using bploc 
tic; p=bploc(p); toc % localize 1st BP 
p=swibra(dir,'bpt1','b',0.01);  p=cont(p,1); % compute 1 point on bif. branch
%% C3: comparison of lss by cont till lam=1; 
p0=loadp('b','pt1'); p0.file.smod=0; p0.plot.pmod=0; 
p=p0; t1=tic; p=cont(p); t1=toc(t1); % using lssbel with lssAMG 
p=p0; p.fuha.lss=@lss; p.fuha.blss=@lss; % switching back to \ 
t2=tic; p=cont(p); t2=toc(t2); fprintf('t1=%g, t2=%g\n',t1,t2);
%% C4: plotting 
figure(3); clf; plotbra('b','lab',12); plotsol('b','pt12',1,1,2);
%% times 
nu=[8000 27000 64000 125000]; lss=[7 73 330 1430]; ilu=[2.6 11 31 72]; 
figure(1); plot(nu,10*ilu,'-*', nu,lss,'-*'); axis tight; legend('10*AMG','lss');
xlabel('n_u'); ylabel('time (seconds)'); 