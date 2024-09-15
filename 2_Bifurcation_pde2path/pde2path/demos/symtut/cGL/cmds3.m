%% hom. NEUMANN BC: continue zero solution
p=loadp('zero','pt0','nbc/zero'); p=setlam(p,-0.1); p.nc.dsmax=0.1; p=cont(p,10);
%% switch to front-mode for one step
p=swibra('nbc/zero','bpt2','nbc/ini',0.01); p.sw.bifcheck=0; p=cont(p,1);
%% add phase conditions and continue
p=loadp('nbc/ini','pt1','nbc/stand'); p=resetc(p); 
p.nc.ilam=[1;2;3]; p.nc.nq=2; p.fuha.qf=@qf_sym;
p.sw.qjac=1; p.fuha.qfder=@qfder_sym; 
p.sw.bprint=[2;3]; clf(2); p.nc.dsmax=0.5; p=cont(p,20); 
%% continue in family given by rotation symmetry
p=loadp('nbc/stand','pt20','nbc/rot'); p=resetc(p); clf(2);
p.nc.ilam=[8;2;3]; p=cont(p,20);  % L2-norm is of the first component only
%% increase domain size and refine mesh
p=loadp('nbc/stand','pt20','nbc/dom'); p=resetc(p); clf(2);
p.nc.ilam=[9;2;3]; p.plot.bpcmp=9;
p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.lammin=0.5; p=cont(p,10);
p=meshada(p,'ngen',3,'sig',1e-4); 
p.nc.dsmax=0.1; p.nc.lammin=0.25; p=cont(p,10);
%% nonzero frequency/speed through gamma
p=loadp('nbc/dom','pt10','nbc/move-'); p=resetc(p); clf(2);
p.nc.ilam=[7;2;3]; p.branch=[bradat(p); p.fuha.outfu(p,p.u)]; 
p.nc.lammin=-1; p.plot.bpcmp=2; p=cont(p,22);
p=loadp('nbc/dom','pt10','nbc/move+'); p=resetc(p); p.sol.ds=-p.sol.ds; 
p.nc.ilam=[7;2;3]; p.branch=[bradat(p); p.fuha.outfu(p,p.u)]; 
p.nc.lammin=-1; p.plot.bpcmp=2; p=cont(p,22);
