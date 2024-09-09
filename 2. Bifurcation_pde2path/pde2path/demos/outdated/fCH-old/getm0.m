function m0=getm0(p) 
m0=triint(p.u(1:p.np),p.mesh.p,p.mesh.t)/p.vol;