function p=qxyron(p) % switch on x-phase condition 
p.nc.nq=3; p.fuha.qf=@qxyr; p.fuha.qfder=@qxyrder; p.nc.ilam=[1 4 5 6];