function out=chembra(p,u)
% outfu for chemotax, put ||u-1||_L^1 first on branch 
n=triint(abs(1-u(1:p.nu)),p.mesh.p,p.mesh.t)/p.vol;
out=[n; max(u(1:p.nu))];
