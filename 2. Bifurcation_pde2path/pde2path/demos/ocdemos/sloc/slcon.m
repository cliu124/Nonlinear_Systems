function q=slcon(p,u) % compute control from states/costates 
q=-1./u(p.np+1:p.nu); 