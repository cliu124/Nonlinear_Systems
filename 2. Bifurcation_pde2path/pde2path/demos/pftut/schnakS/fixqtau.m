function p=fixqtau(p)
stau=size(p.tau,1); p.tau=[p.tau; zeros(p.nu+p.nc.nq+1-stau,1)]; 