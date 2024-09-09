function qu=qfder(p,u)  
qu = (p.mat.Dphi*p.u(1:p.nu))';
end
