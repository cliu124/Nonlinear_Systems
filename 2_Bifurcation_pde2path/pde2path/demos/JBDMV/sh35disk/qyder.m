function qu=qyder(p,u) % derivative of transl.phase condition in y 
qy=(p.mat.Dy(1:p.np,1:p.np)*p.u(1:p.np))';  qu=[qy, 0*qy]; 