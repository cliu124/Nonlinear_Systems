function [c,a,f,b]=gpG(p,u)  % pde for GP after Lashkin08
par=u(p.nu+1:end); lam=par(1); pa=par(2); mu=par(3); c=[1;0;0;1;1;0;0;1]; 
x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
xi=pdeintrp(p.mesh.p,p.mesh.t,x); yi=pdeintrp(p.mesh.p,p.mesh.t,y);
u=p.mat.fill*u(1:p.nu); 
ui=pdeintrp(p.mesh.p,p.mesh.t,u); u=ui(1,:); v=ui(2,:); ua=u.^2+v.^2; 
pot=pa*p.mat.poti;
f1=ua.*u; f2=ua.*v; f=[f1;f2]; a=[-mu+pot;-mu+pot]; % a is diagonal
b=zeros(p.nc.neq^2*2,p.mesh.nt); 
b(3,:)=yi;b(4,:)=-xi;b(5,:)=-yi;b(6,:)=xi; b=lam*b; 
