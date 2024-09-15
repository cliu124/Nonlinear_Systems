function q=qf(p,u) 
uold=p.u(1:p.nu); uox=p.mat.Dphi*uold; %uox=p.mat.Dx*uold; % formally more correct! 
q=uox'*(u(1:p.nu)-uold);