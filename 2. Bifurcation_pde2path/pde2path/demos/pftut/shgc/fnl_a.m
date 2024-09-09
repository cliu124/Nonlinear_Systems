function fa=fnl_a(p,u)  % \pa_a fnl(u,a), where fnl(u,a)=fn(u,<h(u)>) 
ga=p.u(p.nu+3); u1=u(1:p.np); fa=[-ga*u1; 0*u1]; 