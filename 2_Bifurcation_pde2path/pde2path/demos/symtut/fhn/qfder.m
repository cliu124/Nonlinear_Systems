function qu=qfder(p,u)  
qu = (p.mat.Kx*p.u(1:p.nu))';
end
