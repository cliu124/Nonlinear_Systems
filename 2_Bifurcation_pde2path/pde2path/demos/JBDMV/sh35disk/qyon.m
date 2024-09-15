function p=qyon(p) % switch on y-phase condition 
p.nc.nq=1; p.fuha.qf=@qy; p.fuha.qfder=@qyder; p.nc.ilam=[1 5];