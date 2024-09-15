function j=fnljac(p,u)  % \pa_u fnl(u,a), here fnl=-lam*a*u 
h=hfu(p,u); a=p.avvec*h; n=p.np; del=u(p.nu+4); 
j=spdiags(-del*a*ones(n,1),0,n,n); 