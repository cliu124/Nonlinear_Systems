function f=fnl(p,u)  % f(u,<h(u)>) for global coupling 0=G(u)+fnl(u,<h(u)>)
h=hfu(p,u); a=p.avvec*h; ga=u(p.nu+3); u1=u(1:p.np); f=[-ga*a*u1; 0*u1]; 