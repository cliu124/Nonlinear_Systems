function fa=fnl_a(p,u)  % \pa_a fnl(u,a), where fnl(u,a)=fn(u,<h(u)>) 
lam=u(p.nu+1); fa=-lam*u(1:p.np); %fa=ones(p.np,1); 