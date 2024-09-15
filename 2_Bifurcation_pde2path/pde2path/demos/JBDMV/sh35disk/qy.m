function q=qy(p,u) % phase condition for transl.invariance in y 
n=p.np; u0y=p.mat.Dy(1:n,1:n)*p.u(1:n); q=u0y'*u(1:n);