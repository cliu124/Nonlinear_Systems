function [cj,aj,bj]=nlbjac(p,u) 
% jacobian for GP
cj=[1;0;0;1;1;0;0;1]; bj=zeros(8,1); kstar=u(p.nu+2:p.nu+3);
bj(5)=-2*kstar(1);bj(3)=2*kstar(1);bj(6)=-2*kstar(2);bj(4)=2*kstar(2);
lam=u(p.nu+1);ksq=kstar(1)^2+kstar(2)^2;
ufull=p.mat.fill*u(1:p.nu);
ui=pdeintrp(p.mesh.p,p.mesh.t,ufull);
ure=ui(1,:); uim=ui(2,:); 
ua=ure.^2+uim.^2; 
uau=2*ure; uav=2*uim; 
f1u=-p.sig*(ua+uau.*ure)-(ksq+p.mat.poti-lam); f1v=-p.sig*uav.*ure; 
f2u=-p.sig*uau.*uim; f2v=-p.sig*(ua+uav.*uim)-(ksq+p.mat.poti-lam); 
aj=-[f1u; f2u; f1v; f2v]; 
