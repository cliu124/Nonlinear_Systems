function p=qroton(p) % switch on rotational phase-condition 
p.nc.nq=1; p.fuha.qf=@qrot; p.fuha.qfder=@qrotder; p.nc.ilam=[1 6];