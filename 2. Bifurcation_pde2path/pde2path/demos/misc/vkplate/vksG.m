function r=vksG(p,u) 
% von karman-plate as a 10 component system including some regularization
% i.e., aux variables like u5 defined as u5=(1-del Lap)^(-1)\pa_x^2 u1 
lam=u(p.nu+1); 
u5=u(4*p.np+1:5*p.np); u6=u(5*p.np+1:6*p.np); 
u7=u(6*p.np+1:7*p.np); u8=u(7*p.np+1:8*p.np); 
u9=u(8*p.np+1:9*p.np); u10=u(9*p.np+1:10*p.np); 
f1=-(u5.*u9-2*u7.*u10+u6.*u8); 
f2=u5.*u6-u7.^2;  zv=zeros(p.np,1);
f=[zv;f1;zv;f2;zv;zv;zv;zv;zv;zv];
r=(p.mat.K+lam*p.mat.K2)*u(1:p.nu)-p.mat.M*f; 
