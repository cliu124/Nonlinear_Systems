function p=stan2rot(p,tol,ds,mm) % switch on 2 rot.PCs (convenience function) 
p=rotax(p,mm);
p.fuha.qf=@qAV2rot; p.fuha.qfder=@qAV2rotder; p.sw.para=2; p.nc.bisecmax=10;
p.nc.ilam=[6 2 3 7 8 9 10 11];  p.nc.nq=length(p.nc.ilam)-1; p.sw.bifcheck=2; 
p.nc.neig=5; p.nc.eigref=-30; p.nc.tol=tol; p.sol.ds=ds;