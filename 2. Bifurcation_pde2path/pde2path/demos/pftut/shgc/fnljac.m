function j=fnljac(p,u)  % \pa_u fnl(u,<h(u)>) 
h=hfu(p,u); a=p.avvec*h; ga=p.u(p.nu+3); n=p.np; 
j=spdiags([-ga*a*ones(n,1); zeros(n,1)],0,2*n,2*n); 