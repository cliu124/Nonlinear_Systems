function p=qxyon(p) % switch on x-phase condition 
p.nc.nq=2; p.fuha.qf=@qxy; p.fuha.qfder=@qxyder; p.nc.ilam=[1 4 5];