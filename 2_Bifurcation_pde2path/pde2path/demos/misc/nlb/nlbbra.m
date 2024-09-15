function out=nlbbra(p,u)
% output to bifurcation diagram function for nlb 
% max |u|, ||u||_2
u_full = p.mat.fill*u(1:p.nu);
out=[max(abs(u_full(1:p.np)+1i*u_full(p.np+1:2*p.np)));
    sqrt(triint(u_full(1:p.np).^2+u_full(p.np+1:2*p.np).^2,p.mesh.p,p.mesh.t))];