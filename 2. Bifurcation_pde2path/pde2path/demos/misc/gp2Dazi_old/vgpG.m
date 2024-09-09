function [c,a,f,b]=vgpG(p,u) 
% pde for 2-compo GP after Lashkin-Kivshar09, i.e. 4 compo real system
par=u(p.nu+1:end); lam=par(1); pa=par(2); mu1=par(3); mu2=par(4); 
c=[1;0;0;1;1;0;0;1;1;0;0;1;1;0;0;1];
pot=pa*p.mat.poti; x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
xi=pdeintrp(p.mesh.p,p.mesh.t,x);yi=pdeintrp(p.mesh.p,p.mesh.t,y);
u=p.mat.fill*u(1:p.nu); ui=pdeintrp(p.mesh.p,p.mesh.t,u);
u1=ui(1,:); v1=ui(2,:);u2=ui(3,:); v2=ui(4,:); 
a1=u1.^2+v1.^2; a2=u2.^2+v2.^2; 
f1=(a1+0.5*a2).*u1;f2=(a1+0.5*a2).*v1;
f3=(a2+0.5*a1).*u2;f4=(a2+0.5*a1).*v2;
f=[f1;f2;f3;f4]; 
a=[-mu1+pot;-mu1+pot;-mu2+pot;-mu2+pot]; %a diagonal
b=zeros(p.nc.neq^2*2,p.mesh.nt); 
b(3,:)=yi;b(4,:)=-xi;b(9,:)=-yi;b(10,:)=xi;
b(23,:)=yi;b(24,:)=-xi;b(29,:)=-yi;b(30,:)=xi;
b=lam*b; 
