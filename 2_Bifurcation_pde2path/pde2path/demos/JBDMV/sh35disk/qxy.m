function q=qxy(p,u) % x and y phase conditions 
np=p.np; u=u(1:np); uold=p.u(1:np); u0x=p.mat.Dx(1:np,1:np)*uold; 
u0y=p.mat.Dy(1:np,1:np)*uold; q=[u0x'*u; u0y'*u]; 
