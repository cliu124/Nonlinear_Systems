function q=schnak_qf(p,u) 
uox=p.mat.Kx*p.u(1:p.nu); q=uox'*u(1:p.nu);

