function qu=qxyder(p,u) % derivative of all x-y phase conditions 
np=p.np; uold=p.u(1:np); u0x=p.mat.Dx(1:np,1:np)*uold; u0y=p.mat.Dy(1:np,1:np)*uold; 
qu=[u0x' 0*u0x';u0y' 0*u0y']; 