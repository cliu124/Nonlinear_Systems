function out=fchbra(p,u,lam)
% output to bifurcation diagram function for fch 
% mass, max u_1, min u_1, |u_1|_L^2
m=triint(u(1:p.np),p.mesh.p,p.mesh.t)/p.vol; 
%try isreal(p.m0); catch p.m0=m; end
out=[m; u(p.nu+1); u(p.nu+2); u(p.nu+3); % here eta1,gamma,eps
    max(u(1:p.np)); min(u(1:p.np));...
     sqrt(triint(u(1:p.np).^2,p.mesh.p,p.mesh.t))];
end