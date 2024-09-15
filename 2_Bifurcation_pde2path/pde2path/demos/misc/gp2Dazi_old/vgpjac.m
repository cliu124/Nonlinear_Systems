function [cj,aj,bj]=vgpjac(p,u) 
% jacobian for vector GP 
par=u(p.nu+1:end); lam=par(1); pa=par(2); mu1=par(3); mu2=par(4); 
cj=[1;0;0;1;1;0;0;1;1;0;0;1;1;0;0;1]; 
pot=pa*p.mat.poti; x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
xi=pdeintrp(p.mesh.p,p.mesh.t,x);yi=pdeintrp(p.mesh.p,p.mesh.t,y);
u=p.mat.fill*u(1:p.nu); ui=pdeintrp(p.mesh.p,p.mesh.t,u); 
u1=ui(1,:); v1=ui(2,:);u2=ui(3,:); v2=ui(4,:); 
a1=u1.^2+v1.^2; a2=u2.^2+v2.^2; 
a1u1=2*u1; a1v1=2*v1; a2u2=2*u2; a2v2=2*v2; 
f1u1=mu1-pot+a1u1.*u1+a1+0.5*a2; 
f1v1=a1v1.*u1; f1u2=0.5*a2u2.*u1; f1v2=0.5*a2v2.*u1; 
f2v1=mu1-pot+a1v1.*v1+a1+0.5*a2; 
f2u1=a1u1.*v1; f2u2=0.5*a2u2.*v1; f2v2=0.5*a2v2.*v1; 
f3u2=mu2-pot+a2u2.*u2+a2+0.5*a1;
f3u1=0.5*a1u1.*u2; f3v1=0.5*a1v1.*u2; f3v2=a2v2.*u2; 
f4v2=mu2-pot+a2v2.*v2+a2+0.5*a1;
f4u1=0.5*a1u1.*v2; f4v1=0.5*a1v1.*v2; f4u2=a2u2.*v2;
fu=[f1u1; f2u1; f3u1; f4u1; f1v1; f2v1; f3v1; f4v1; ...
    f1u2; f2u2; f3u2; f4u2; f1v2; f2v2; f3v2; f4v2]; 
bj=zeros(p.nc.neq^2*2,p.mesh.nt); 
bj(3,:)=yi;bj(4,:)=-xi;
bj(9,:)=-yi;bj(10,:)=xi;
bj(23,:)=yi;bj(24,:)=-xi;
bj(29,:)=-yi;bj(30,:)=xi;
bj=lam*bj; aj=-fu; 

