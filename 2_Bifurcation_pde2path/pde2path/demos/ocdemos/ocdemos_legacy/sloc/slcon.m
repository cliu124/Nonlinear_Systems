function k=slcon(p,u) % extract control from states/costates 
k=-1./u(p.np+1:p.nu); 