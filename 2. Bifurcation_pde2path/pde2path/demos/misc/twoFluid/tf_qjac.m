function qu=tf_qjac(p,u)  
% define auxiliary function: here phase condition
uold=p.mat.fill*p.u(1:p.nu);
[~,uoy]=pdegrad(p.mesh.p,p.mesh.t,uold); % y-derivative
ugn=pdeprtni(p.mesh.p,p.mesh.t,uoy); % nodal 
ug=reshape(ugn,p.nc.neq*p.np,1); qu=(p.mat.drop*ug)';
end