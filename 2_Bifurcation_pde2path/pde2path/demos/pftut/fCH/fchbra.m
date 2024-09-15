function out=fchbra(p,u) % output to bifurcation diagram function for fch 
% mass, max u_1, min u_1, |u_1|_L^2
u1=u(1:p.np); m=sum(p.mat.Ms*u1)/p.vol; 
%try isreal(p.m0); catch p.m0=m; end
out=[m; u(p.nu+1:end); % [eta1; ga; eps; mi; eta2]; 
    max(u(1:p.np)); min(u(1:p.np));...
    sqrt(sum(p.mat.Ms*(u1.^2)))];
end