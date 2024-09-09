function [K,M]=LBtor6(p,R,rho)  
% LBtor6: Laplace-Beltrami-Op for torus based on 6-node triangles
% hence assem6g instead assema, otherwise completely as LBtor. 
po=getpte(p); th=po(2,:)'; n=p.np; [Kphi,M]=assem6g(p,[1 0; 0 0],1); 
[Kth,~]=assem6g(p,[0*th' 0*th'; 0*th' R+rho*cos(th)'],1); 
M=filltrafo(p,M);  dd=1./(R+rho*cos(th)); % build LBO and tranform to pBC: 
K=filltrafo(p,spdiags(dd/rho^2,0,n,n)*Kth+spdiags(dd.^2,0,n,n)*Kphi); 