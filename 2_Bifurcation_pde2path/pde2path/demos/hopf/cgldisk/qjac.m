function qu=qjac(p,u) % phase condition jac for only translation
qu=p.u0x(1:p.nu)'; 
end

