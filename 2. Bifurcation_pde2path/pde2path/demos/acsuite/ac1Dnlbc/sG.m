function r=sG(p,u)  % rhs for ac1Dnlbc
f=nodalf(p,u); par=u(p.nu+1:end); d=par(4); al=par(5); beta=par(6); u=u(1:p.nu);
r=par(1)*p.mat.K*u-p.mat.M*f...   % bulk part of PDE 
  +p.mat.Q1*u-beta*p.mat.G1 ... % left BCs, i.e., Q1,G1=0 at right 
  +p.nc.sf*(p.mat.Q2*(u+al*u.^2)-d*p.mat.G2); % right BCs