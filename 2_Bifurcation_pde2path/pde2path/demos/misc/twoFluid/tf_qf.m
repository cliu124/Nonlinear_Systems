function q=tf_qf(p,u)  
% define auxiliary function: here phase condition
uold=p.mat.fill*p.u(1:p.nu);
[~,uoy]=pdegrad(p.mesh.p,p.mesh.t,uold); 
uoyn=pdeprtni(p.mesh.p,p.mesh.t,uoy); % back to nodes!
uoyn=reshape(uoyn,p.np*p.nc.neq,1); 
q= (p.mat.drop*uoyn)'*u(1:p.nu);
end
