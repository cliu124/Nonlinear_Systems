function [K,M]=LBcone6(p,a)  
% LBcone: LB for cone of height 1 and 'circular' base of radius a 
fem=p.pdeo.fem; gr=p.pdeo.grid; po=getpte(p); 
x=po(1,:); y=po(2,:); r2=x.^2+y.^2; % r2=max(rmin,x.^2+y.^2);  % extract coordinates 
c1=a^2+y.^2./r2; c2=-x.*y./r2; c3=c2; c4=a^2+x.^2./r2; % diffusion coeffients 
cc=[[c1 c2]; [c3 c4]]; % stored as a "2 x (2*np) matrix", see grid2D.aCoefficientsMpt
[K,M,~]=fem.assema(gr,cc,1,0); K=K/(a^2*(1+a^2)); 