function r=sG(p,u)  % rhs for ql-AC with pBCs 
par=u(p.nu+1:end); c0=par(1); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
u=u(1:p.nu); uf=p.mat.fill*u; ut=(p.mat.p2c*uf)'; % interpol. to elem. centers
c=c0+del*ut+epsi*ut.^2; f=lam*p.xft.*ut+ut.^3-ga*ut.^5; % diff. tensor and f
[K,~,F]=p.pdeo.fem.assema(p.pdeo.grid,c,0,f); % assemble K and F (M not used) 
K=filltrafo(p,K); F=p.mat.fill'*F; r=K*u-F; % reduce K,F, then compute resi 