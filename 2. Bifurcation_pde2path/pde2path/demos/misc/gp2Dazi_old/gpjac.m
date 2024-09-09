function [cj,aj,bj]=gpjac(p,u) 
% jacobian for GP
par=u(p.nu+1:end); lam=par(1); pa=par(2); mu=par(3); cj=[1;0;0;1;1;0;0;1]; 
pot=pa*p.mat.poti; x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)';
xi=pdeintrp(p.mesh.p,p.mesh.t,x); yi=pdeintrp(p.mesh.p,p.mesh.t,y);
u=p.mat.fill*u(1:p.nu); 
ui=pdeintrp(p.mesh.p,p.mesh.t,u); u=ui(1,:); v=ui(2,:); ua=u.^2+v.^2; 
g=ua; gu=2*u; gv=2*v; 
f1u=mu-pot+gu.*u+g; f1v=gv.*u; 
f2u=gu.*v; f2v=mu-pot+gv.*v+g; 
fu=[f1u; f2u; f1v; f2v]; aj=-fu; 
bj=zeros(p.nc.neq^2*2,p.mesh.nt); 
bj(3,:)=yi;bj(4,:)=-xi;bj(5,:)=-yi;bj(6,:)=xi;bj=lam*bj; 
