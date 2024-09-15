function [c,a,f,b]=ac_G(p,u)  
% nodal version of PDE residual for AC model
    par=u(p.nu+1:end); u=u(1:p.nu); 
    u=pdeintrp(p.mesh.p,p.mesh.t,u);
    f=par(1)*u+u.^3-par(3)*u.^5; 
    a=0;b=0;
    c=par(2);
end