function [c,a,f,b]=schnakf(p,u) 
% Schnakenberg: u_t=Del u-u+u^2v,  v_t=d Del v+lam-u^2 v
lam=u(p.nu+1); u=u(1:p.nu); [po,tri]=getpte(p); u=pdeintrp(po,tri,u); 
d1=1; d2=60; c=[d1;0;0;d1;d2;0;0;d2]; a=0; 
f1=-u(1,:)+u(1,:).^2.*u(2,:); f2=lam-u(1,:).^2.*u(2,:); f=[f1;f2]; b=0;
