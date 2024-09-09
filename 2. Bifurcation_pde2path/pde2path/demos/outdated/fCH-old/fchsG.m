function r=fchsG(p,u)  
% pearling for fCH, nodal version 
% -eps^2 Del u+W'(u)+v=0
% -eps^2 Del v+W''(u)v-eps*eta1*v-eps*etad*W'(u)-eps*ga=0 
f=nf(p,u); par=u(p.nu+1:end); eps=par(3);
r=eps^2*p.mat.K*u(1:p.nu)-p.mat.M*f; 
end

function f=nf(p,u)
par=u(p.nu+1:end); u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
eta1=par(1); ga=par(2); eps=par(3); eta2=par(5); etad=eta1-eta2; 
[w,wp,wpp,wppp]=p.fuha.wfu(u1,p); 
f1=-wp-u2; 
f2=-wpp.*u2+eps*eta1*u2+eps*etad*wp-eps*ga; 
f=[f1;f2]; 
end