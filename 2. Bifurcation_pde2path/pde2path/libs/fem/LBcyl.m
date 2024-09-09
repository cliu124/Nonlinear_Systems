function K=LBcyl(p,u,R) 
% LBcyl: Laplace-Beltrami for cylinder of radius R, coord (phi,z)
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[Kphi,M,~]=fem.assema(gr,[ones(n,1);0*ones(n,1)],1,1); 
[Kz,~,~]=fem.assema(gr,[0*ones(n,1); ones(n,1)],0,1); 
[Kphi,~,~]=fem.assema(gr,[1 0; 0 0],1,1); 
[Kth,~,~]=fem.assema(gr,[0 0; 0 1],0,1); 
K=filltrafo(p,spdiags(dd/R^2,0,n,n)*Kth+spdiags(dd.^2,0,n,n)*Kphi); 