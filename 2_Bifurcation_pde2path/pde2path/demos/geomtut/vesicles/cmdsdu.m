%% 3rd from oblate 
p=swibra('2/o','bpt3','2/duo-3',-0.01); p.sw.bifcheck=0; p.nc.tol=1e-3; % no bifcheck, poor tol 
p.nc.Ab=inf; p.nc.delbound=inf; % no area refinement in the first steps
p=cont(p,4); 
%% Full rotational PC's
p=loadp('0/o-1','pt4'); p=stanfullrot(p,1e-8,-0.02); 
p.nc.Ab=3*err; p.nc.delbound=15; p.usrlam=[]; p.nc.eigref=-20; 
p.nc.foldtol=0.05; p.nc.dsmax=0.2; p=cont(p,30);