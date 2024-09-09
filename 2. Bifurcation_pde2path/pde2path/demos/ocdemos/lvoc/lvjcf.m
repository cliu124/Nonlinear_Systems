function [jc,jc1,jc2]=lvjcf(p,u) % current value(s) J 
par=u(p.nu+1:end);  al1=par(9); al2=par(10); 
p1=par(11); p2=par(12); c1=par(13); c2=par(14); 
v1=u(1); v2=u(p.np+1); k=lvcon(p,u); 
jc1=p1*v1^al1*k(1)^(1-al1)-c1*k(1); jc2=p2*v2^al2*k(2)^(1-al2)-c2*k(2); 
jc=jc1+jc2; 