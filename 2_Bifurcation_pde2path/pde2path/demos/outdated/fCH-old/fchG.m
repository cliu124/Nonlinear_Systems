function [c,a,f,b]=fchG(p,ug) 
% pearling for fCH
% -eps^2 Del u+W'(u)+v=0
% -eps^2 Del v+W''(u)v-eps*eta1*v-eps*etad*W'(u)-eps*ga=0 
% ga=u(nq+1) 
global eta; % check dimension of eta here, needed in fchqf
u=ug(1:p.nu); par=ug(p.nu+1:end); 
eta1=par(1); ga=par(2); eps=par(3); eta2=par(5); etad=eta1-eta2; 
try se=size(eta,2); catch; eta=[]; se=0; end
if(se~=size(ug,1)) % eta not yet set, or mesh is refined 
  C=n2triamat(p.mesh.p,p.mesh.t); ta=triar(p.mesh.p,p.mesh.t); 
  eta=ta*C; eta=[eta zeros(1,p.np)]; 
end
u=pdeintrp(p.mesh.p,p.mesh.t,u); 
d1=eps^2; d2=eps^2; c=[d1;0;0;d1;d2;0;0;d2]; 
[w,wp,wpp,wppp]=p.fuha.wfu(u(1,:),p); 
a=[0;0;1;-eps*eta1];
f1=-wp; f2=-wpp.*u(2,:)+eps*etad*wp-eps*ga; 
f=[f1;f2]; b=0;
