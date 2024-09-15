function r=sG(p,u)
    u1=u(1:p.np); % identify solution component 1
    u2=u(p.np+1:2*p.np); % identify solution component 2
    par=u(p.nu+1:end); % identify parameters
    f1=-u1+u1.^2.*u2+par(2)*(u1-u2.^(-1)).^2; % non-linearity for u1 
    f2=par(1)-u1.^2.*u2-par(2)*(u1-u2.^(-1)).^2; % non-linearity for u2
    f=[f1;f2];
    r=kron([[1,0];[0,par(3)]],p.mat.K)*u(1:p.nu)-p.mat.M*f; % calculation of the
                                                            % residual
end


