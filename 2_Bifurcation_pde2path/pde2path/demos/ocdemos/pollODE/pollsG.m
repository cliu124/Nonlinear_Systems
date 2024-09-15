function r=pollsG(p,u) % rhs for pollution problem 
f=nodalf(p,u); r=-p.mat.M*f; % the residual 