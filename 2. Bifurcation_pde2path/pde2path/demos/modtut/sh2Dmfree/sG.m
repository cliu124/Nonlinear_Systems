function r=sG(p,u)  % SH 2D via dct with preassemble dct-matrix F and mult. mu
global p2pglob; F=p2pglob.F; 
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(1); c2=par(2); c3=par(3); % split 
u=F'*uf; ff=lam*u+c2*u.^2+c3*u.^3;  f=F*ff; % "nonlinearity"  
r=p2pglob.mu.*uf-f;     % residual 