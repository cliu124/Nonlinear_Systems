function q=qf(p,u) 
uold=p.u0; uox=p.u0x;  q=uox'*(u(1:p.nu)-uold);
end