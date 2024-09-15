function qu=acfront_qfder(p,u)  
uold=p.mat.fill*p.u(1:p.nu);
[~,uoy]=pdegrad(p.mesh.p,p.mesh.t,uold); 
% u-derivative of auxiliary function: here phase condition
uyn=pdeprtni(p.mesh.p,p.mesh.t,uoy); % nodal
qu=reshape(uyn,p.np,1)';
end