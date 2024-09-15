function [K,M]=LBcone2(p,a,e)  
% LBcone2: LB for cone of height a and elliptic (e,1) plane base  
fem=p.pdeo.fem; gr=p.pdeo.grid; po=getpte(p); e2=e^2; a2=a^2; 
x=po(1,:); y=po(2,:); r2=x.^2/e2+y.^2; % extract coordinates 
d=1+a^2*(x.^2/(e2^2)+y.^2)./r2; sdi=1./sqrt(d);
c1=sdi.*(1+a2*y.^2./r2); 
c2=-a2*sdi.*x.*y./(e2*r2); 
c3=c2; c4=sdi.*(1+a2*x.^2/(e2^2*r2)); % diffusion coeffients 


cc=[[c1 c2]; [c3 c4]]; % stored as a "2 x (2*np) matrix", see grid2D.aCoefficientsMpt
[K,M,~]=fem.assema(gr,cc,1,0); SD=spdiags(sdi',0,p.np,p.np); K=SD*K; 