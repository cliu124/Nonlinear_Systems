function qu=qrotder(p,u) % derivative of rotational phase-cond 
n=p.np; qr=(p.mat.Krot(1:n,1:n)*p.u(1:n))'; qu=[qr, 0*qr]; 