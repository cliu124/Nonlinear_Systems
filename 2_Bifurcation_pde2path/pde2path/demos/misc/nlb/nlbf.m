function [c,a,f,b]=nlbf(p,u) 
c=[1;0;0;1;1;0;0;1]; b=zeros(8,1);
kstar=u(p.nu+2:p.nu+3);
lam=u(p.nu+1);
b(5)=-2*kstar(1);b(3)=2*kstar(1);b(6)=-2*kstar(2); b(4)=2*kstar(2);
ufull=p.mat.fill*u(1:p.nu); 
ui=pdeintrp(p.mesh.p,p.mesh.t,ufull); ure=ui(1,:); uim=ui(2,:); ua=ure.^2+uim.^2; 
f1=-p.sig*ua.*ure; f2=-p.sig*ua.*uim; f=[f1;f2]; 
ksq = kstar(1)^2+kstar(2)^2;
a=[ksq+p.mat.poti-lam; ksq+p.mat.poti-lam]; % a is diagonal