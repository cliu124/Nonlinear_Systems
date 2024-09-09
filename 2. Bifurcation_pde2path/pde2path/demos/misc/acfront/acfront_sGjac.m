function Gu=acfront_sGjac(p,u)  
    par=u(p.nu+1:end);
    fu=par(1)*par(2)-3*par(1)*u.^2 + 2*par(1)*(1-par(2))*u; %ac jacobian nodal version
    Fu=spdiags(fu,0,p.nu,p.nu);
    Gu=p.mat.K-par(3)*p.mat.Kadv-p.mat.M*Fu; 
end
