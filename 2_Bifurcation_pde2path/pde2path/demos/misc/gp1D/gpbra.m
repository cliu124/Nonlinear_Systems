function out=gpbra(p,u) % GPbra 
nu=p.nu; np=p.np; par=u(nu+1:end); u1=u(1:nu); 
N=sum(p.mat.M*(u1.^2)); H=ham(p,u); 
out=[par; max(abs(u1(1:np))); N; H]; 