function z=zfu(p)
np=p.np; z=exp(p.u(1:np));
z=p.u(np+1:2*np);