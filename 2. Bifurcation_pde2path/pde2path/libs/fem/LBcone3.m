function [K,M]=LBcone3(p,h,a)  
% LBcone3: LB cone, height h and ellipticity a, parametr. over unit disk
fem=p.pdeo.fem; gr=p.pdeo.grid; po=getpte(p); h2=h^2; a2=a^2; 
x=po(1,:); y=po(2,:); r2=x.^2+y.^2; % extract coordinates 
d=1+a2*(x.^2+a2*y.^2)./r2; sd=sqrt(d); sdi=1./sd;
c1=sdi.*(1+h2*y.^2./r2); 
c2=-h2*sdi.*x.*y./r2; 
c3=c2; c4=sdi.*(a2+h2*x.^2./r2); % diffusion coeffients 


cc=[[c1 c2]; [c3 c4]]; % stored as a "2 x (2*np) matrix", see grid2D.aCoefficientsMpt
[K,M,~]=fem.assema(gr,cc,1,0); SD=spdiags(sd',0,p.np,p.np); %SDI=spdiags(sd',0,p.np,p.np); 
%K=SD*K;  
M=SD*M; 