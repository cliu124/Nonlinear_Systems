function r=sG(p,u)  % AC 1D with dct-diff.matrix, u in x-space 
n=p.nu; par=u(n+1:end); u=u(1:n); % split u into par and pde parts 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % nonlin. 
r=par(1)*p.mat.L*u-f;    % bulk part of PDE 
  