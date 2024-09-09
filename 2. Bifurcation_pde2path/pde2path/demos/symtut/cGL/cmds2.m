%% fold-continuation 
p=spcontini('per/stand','fpt2',6,'per/fc1'); % init fold cont, par 6 new prim. par
p.plot.bpcmp=1; figure(2); clf; p.sol.ds=-0.01; % use this new param.for plotting
p.sw.spjac=1; p.fuha.spjac=@spjac; % spectral jac
p.sw.spqjac=1; p.fuha.spqjac=@spqjac; % and jac for constraints
tic; p=cont(p); toc
%% switch back to regular cont
p=spcontexit('per/fc1','pt20','per/stand2a'); p.plot.bpcmp=0; 
clf(2); p=cont(p,1); p=cont(p,20); % cont in one direction, 1 init.step for saving
% continue in other direction 
p=loadp('per/stand2a','pt1','per/stand2b'); p.sol.ds=-p.sol.ds; p=cont(p,20); 
%% continue standing wave in family given by rotation symmetry
p=loadp('per/stand','pt20','per/rot'); p=resetc(p); p.file.smod=1; 
p.nc.ilam=[8;2;3]; p.nc.dsmax=0.3; p.nc.lammax=10; p=cont(p,30);
%% continue WT in family given by rot=transl symmetry (shape doesn't change)
p=loadp('per/wtstand','pt20','per/trans'); p=resetc(p); clf(2);
p.nc.ilam=[8;2]; p.nc.lammax=10; p=cont(p,20);
%% continue WT to nonzero speed/rotation (shape doesn't change)
p=loadp('per/wtstand','pt20','per/wtmove'); p=resetc(p); clf(2);
p.nc.ilam=[4;2]; p.sol.ds=-p.sol.ds; p.plot.bpcmp=2; 
p.nc.tol=1e-4; p=cont(p,2); p.nc.tol=1e-8; p=cont(p,20); 
%% break rotation symmetry through gamma for a wavetrain
p=loadp('per/wtstand','pt20','per/wtasym'); p=resetc(p); clf(2);
p.sw.foldcheck=0;  p.nc.ilam=[7;2]; p.plot.bpcmp=0; 
p.nc.tol=1e-5; p=cont(p,2); p.nc.tol=1e-8; p=cont(p,20); 