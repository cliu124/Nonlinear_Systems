function q=qf(p,u) 
uold=p.u(1:p.nu); uox = p.mat.Kx*uold;
q=uox'*(u(1:p.nu)-uold);
end