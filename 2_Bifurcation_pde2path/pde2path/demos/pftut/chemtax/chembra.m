function out=chembra(p,u)
% outfu for chemotax, put ||u-1||_L^1 first on branch 
n=sum(abs(p.mat.M(1:p.np,1:p.np)*(1-u(1:p.np))))/p.vol; 
out=[n; max(u(1:p.nu))];
