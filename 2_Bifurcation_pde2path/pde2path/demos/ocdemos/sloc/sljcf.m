function jc=sljcf(p,u) % current value J 
cp=u(p.nu+3); v=u(1:p.np); q=-1./u(p.np+1:p.nu); jc=log(q)-cp*v.^2;