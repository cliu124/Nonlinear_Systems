function p=qxon(p) % switch on x-phase condition 
p.nc.nq=1; p.fuha.qf=@qx; p.fuha.qfder=@qxder; p.nc.ilam=[1 4];