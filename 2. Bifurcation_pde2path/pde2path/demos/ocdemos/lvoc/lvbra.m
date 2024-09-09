function out=lvbra(p,u) % branch-output for lvoc 
par=u(p.nu+1:end); p1=par(11); p2=par(12); ga1=par(6); ga2=par(7);
v1=u(1); v2=u(p.np+1); l1=u(2*p.np+1); l2=u(3*p.np+1); 
au1=p1-ga1*l1; au2=p2-ga2*l2; [h,hp]=hfu(p,p.u);
[jc,jc1,jc2]=lvjcf(p,p.u);   k=lvcon(p,p.u); 
out=[jc; jc1; jc2; v1; v2; l1; l2; k(1); k(2); h(1); h(2); au1; au2; 
    max(u(1:p.np)); min(u(1:p.np))];