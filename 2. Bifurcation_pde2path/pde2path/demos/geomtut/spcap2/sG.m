function r=sG(p,u)  % PDE rhs 
par=u(p.nu+1:end); H0=par(1); u=u(1:p.nu); al=par(3); % split in par and PDE u 
N0=getN(p,p.X); X=p.X+u.*N0; H=getmeancurv(p,X); % mean curv.based on FEM 
r=-2*H+p.mat.M*(H0*ones(p.nu,1)); r(p.idx)=u(p.idx)-al; % residual, and DBCs 