function [K,M,Ms]=LBSd(p,a,c)  
% LBSd: Laplace-Beltrami (and mass M) for spheroid with semi-axis a,c
fem=p.pdeo.fem; gr=p.pdeo.grid; po=getpte(p); n=p.np; 
y=po(2,:)'; zv=zeros(n,1); 
g11=a^2*cos(y).^2; g22=a^2*sin(y).^2+c^2*cos(y).^2; 
g=g11.*g22; gs=sqrt(g); 
gsurf=sqrt(a^2*cos(y).^2.*(c^2*cos(y).^2+a^2*sin(y).^2)); % surface element 
[Kphi,M,~]=fem.assema(gr,[gs./g11; zv],1,1);  
[Kth,Ms,~]=fem.assema(gr,[zv; gs./g22],gsurf,1); 
K=filltrafo(p,spdiags(1./gs,0,n,n)*(Kphi+Kth)); 
M=filltrafo(p,M); Ms=filltrafo(p,Ms); 


