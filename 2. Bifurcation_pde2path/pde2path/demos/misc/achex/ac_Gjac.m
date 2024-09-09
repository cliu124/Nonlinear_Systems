function [cj,aj,bj]=ac_Gjac(p,u)  
% Nodal version of Jacobian of PDE for AC model.
    par=u(p.nu+1:end); u=u(1:p.nu);
    u=pdeintrp(p.mesh.p,p.mesh.t,u);
    fu=par(1)+3*u.^2-par(3)*5*u.^4; 
    cj=par(2); aj=-fu;
    bj=0;
end
