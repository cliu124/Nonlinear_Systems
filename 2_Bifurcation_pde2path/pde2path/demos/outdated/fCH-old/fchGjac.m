function [cj,afj,bj]=fchGjac(p,ug) 
% Jacobian for fCH, ignores mass-cons-pen; 
u=ug(1:p.nu); par=ug(p.nu+1:end); 
eta1=par(1); ga=par(2); eps=par(3); eta2=par(5); etad=eta1-eta2; 
u=pdeintrp(p.mesh.p,p.mesh.t,u);
d1=eps^2; d2=eps^2; cj=[d1;0;0;d1;d2;0;0;d2]; bj=0;
[w,wp,wpp,wppp]=p.fuha.wfu(u(1,:),p); ov=ones(1,p.mesh.nt); 
f1u=-wpp; f1v=-ov; 
f2u=-wppp.*u(2,:)+eps*etad*wpp; f2v=-wpp+eps*eta1*ov; 
fu=[f1u; f2u; f1v; f2v]; afj=-fu; 