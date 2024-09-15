%% c1: fold-continuation 
p=spcontini('b1','fpt1',3,'b1f');   % init fold cont with par 3 new prim. par
p.plot.bpcmp=p.nc.ilam(2); figure(2); clf; % use this new param. for plotting
p.sol.ds=-0.01;                   % new stepsize in new primary parameter
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
[Ja, Jn]=spjaccheck(p); pause % check impl. of spjac (comment out when fine) 
tic; p=cont(p); toc
clf(3); plotbraf('b1f','pt20',3,2,'lab',15); % plot BD for fold-cont
%% c2: switch back to regular continuation from one of the fold points
p=spcontexit('b1f','pt15','b1-a'); p.nc.dsmax=0.2; p.sw.bifcheck=0; 
p.plot.bpcmp=0; p.nc.lammin=-5; p.sol.ds=1e-3; clf(2); 
p=cont(p,1); p=cont(p,35); % cont. in one direction, 1 initial step for saving
%% c3: continue in other direction 
p=loadp('b1-a','pt1','b1-b'); p.nc.dsmax=0.1; p.sol.ds=-p.sol.ds; p=cont(p,50); 
%% c4: plot new branches
figure(3); clf; f=3; c=0; 
plotbra('b1-a','pt20',f,c,'fplab',1,'lsw',0); 
plotbra('b1-b','pt40',f,c,'fplab',1,'lsw',0); 
plotbra('b1','pt20',f,c,'cl','r','fplab',1,'lsw',0); 
%% c5 - branch point continuation, here in c for first BP on trivial branch 
p=bpcontini('tr','bpt1',1,'bpcon'); % init branch point cont. in par(1) (c)
p.sw.bifcheck=0; p.plot.bpcmp=2; % plot lam over c now
p.sol.ds=-0.1; p.fuha.spjac=@spjac; p.usrlam=0.3; clf(2); p=cont(p,10); 
%% c6 - plot lam-pos over c
figure(3); clf; plotbra('bpcon', 'pt10');
%% c7 - exit BP cont 
p=bpcontexit('bpcon','pt7','trb'); % switch to continuation in lambda again
%% c8 - Branch-sw. at continued BP
p=swibra('trb','bpt1','b1b',0.02); % switch to non-trivial branch
p.sw.spcalc=1; p.nc.dsmax=0.15; clf(2); p.usrlam=0; p=cont(p,20); 
%% c9 - plots 
fnr=3; figure(fnr); clf; 
plotbra('b1','pt10',fnr,0,'lab',8);plotbra('b1b','pt20',fnr,0,'cl','r','lab',17); 
plotsol('b1','pt11',1,1,1); pause; plotsol('b1b','pt17',1,1,1);