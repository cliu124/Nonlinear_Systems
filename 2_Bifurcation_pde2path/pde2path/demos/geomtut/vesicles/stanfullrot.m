function p=stanfullrot(p,tol,ds) % switch on full rotational PCs (convenience function) 
p.fuha.qf=@qAVfrot; p.fuha.qfder=@qAVfrotder; p.sw.para=2; p.nc.bisecmax=10;
p.nc.ilam=[6 2 3 7 8 9 10 11 12];  p.nc.nq=length(p.nc.ilam)-1; p.sw.bifcheck=2; 
p.nc.neig=5; p.nc.eigref=0; p.nc.tol=tol; p.sol.ds=ds;