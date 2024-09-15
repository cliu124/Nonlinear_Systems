function f=fnl(p,u)  % f(u,a) for GC 
h=hfu(p,u); a=p.avvec*h; del=u(p.nu+4); f=-del*a*u(1:p.np);  