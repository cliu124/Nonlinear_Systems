function q=acfront_qf(p,u)  
uold=p.mat.fill*p.u(1:p.nu);
[~,uoy]=pdegrad(p.mesh.p,p.mesh.t,uold); 
uoy=pdeprtni(p.mesh.p,p.mesh.t,uoy); % back to nodes!
uoy=reshape(uoy,p.np*p.nc.neq,1); 
% phase condition in periodic case: <uold',u>=0
q= uoy'*(u(1:p.nu)-p.u(1:p.nu));

end
