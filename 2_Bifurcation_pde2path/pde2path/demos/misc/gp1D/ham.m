function H=ham(p,u)
par=u(p.nu+1:end); u=u(1:p.nu); s=par(2); sig=par(3); ga=par(4); 
x=getpte(p); x=x'; V=pot(x,s);
d1=(p.mat.K*u).*u; d2=ga/(sig+1)*u.^(2*sig+2)+V.*(u.^2); 
h1=sum(d1,1); h2=sum(p.mat.M*d2,1); H=h1+h2; 