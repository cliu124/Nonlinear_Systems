function n=l1norm(u)
% L1NORM: l1 norm of u, nodal and scalar 
%
% n=l1norm(u) 
n=triint(abs(u)); 
