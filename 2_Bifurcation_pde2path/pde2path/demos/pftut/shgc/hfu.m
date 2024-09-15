function h=hfu(p,u) % h for GC
u=u(1:p.np); h=[u.^2; zeros(p.np,1)];  
