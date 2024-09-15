function k=cficon(p,u) % extract control from states/costates 
par=u(p.nu+1:end); 
ga1=par(6); ga2=par(7); al1=par(9); al2=par(10); 
p1=par(11); p2=par(12); c1=par(13); c2=par(14); 
v1=u(1); v2=u(p.np+1); l1=u(2*p.np+1); l2=u(3*p.np+1); 
k1=((p1-ga1*l1)*(1-al1)^2/c1)^(1/al1)*v1;
k2=((p2-ga2*l2)*(1-al2)^2/c2)^(1/al2)*v2;
k=[k1 k2]; return
p1-ga1*l1, p2-ga2*l2 %, (p1-ga1*l1)^(1/al1), (p2-ga2*l2)^(1/al2), 1/al1