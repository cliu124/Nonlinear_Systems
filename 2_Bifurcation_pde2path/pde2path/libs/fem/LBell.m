function [K,M]=LBell(p,a,b)  
% LBell: 1D Laplace-Beltrami (and mass M) for ellipse, semiaxis a,b 
fem=p.pdeo.fem; gr=p.pdeo.grid; th=getpte(p); n=p.np; 
dd=a^2*sin(th).^2+b^2*cos(th).^2; dds=1./sqrt(dd); 
[K,M,~]=fem.assema(gr,dds,1,1);  
K=filltrafo(p,spdiags(dds',0,n,n)*K); M=filltrafo(p,M); 