function [K,M]=LBtor(p,R,rho)  
% LBtor: Laplace-Beltrami (and mass M) for torus, radii R and rho 
% (parametric)) domain (phi,th)\in [-\pi,\pi]^2
fem=p.pdeo.fem; gr=p.pdeo.grid; po=getpte(p); th=po(2,:)'; n=p.np; 
[Kphi,M,~]=fem.assema(gr,[1 0; 0 0],1,1); 
[Kth,~,~]=fem.assema(gr,[0*th' 0*th'; 0*th' R+rho*cos(th)'],0,1); 
dd=1./(R+rho*cos(th)); 
K=filltrafo(p,spdiags(dd/rho^2,0,n,n)*Kth+spdiags(dd.^2,0,n,n)*Kphi); 
M=filltrafo(p,M); 