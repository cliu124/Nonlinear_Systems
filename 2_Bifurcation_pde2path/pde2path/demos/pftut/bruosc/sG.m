function r=sG(p,u) % 
f=nodalf(p,u); par=u(p.nu+1:end); % a,b,d_u, d_v,al,om
du=par(3); dv=par(4); 
K=kron([[du,0];[0,dv]],p.mat.K); 
r=K*u(1:p.nu)-p.mat.M*f; 