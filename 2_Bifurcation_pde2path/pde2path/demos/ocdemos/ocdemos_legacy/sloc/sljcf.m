function jc=sljcf(p,u) % current value J 
cp=u(p.nu+3); v=u(1:p.np); k=-1./u(p.np+1:p.nu); jc=log(k)-cp*v.^2;