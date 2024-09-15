function r=pollsG(p,u) % rhs for pollution problem 
f=nodalf(p,u); r=p.mat.K*u(1:p.nu)-p.mat.M*f; % the residual 