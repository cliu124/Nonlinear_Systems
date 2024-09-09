function Gu=sGjac(p,u)  
global p2pglob; % p2pglob.avec, cvec computed here and used in, e.g., gclss
hj=hjac(p,u); p2pglob.avec=sum(p.mat.M*spdiags(hj,0,p.nu,p.nu))./p.Om; 
fnlj=fnljac(p,u); fnla=fnl_a(p,u); p2pglob.cvec=p.mat.M*fnla; 
par=u(p.nu+1:end); u=u(1:p.nu);  % split u into parameters and PDE variables 
K=par(1)*p.mat.K; fu=par(2)+3*u.^2-5*par(3)*u.^4; % K, and local derivatives 
Fu=spdiags(fu,0,p.nu,p.nu); 
% if jfac=0, then lam*<u^j1> in jacobian ignored here and dealt with in (b)gclss.m 
if p.jfac; Fu=Fu+fnla*p2pglob.avec; end % for the brute force impl. of GC terms 
Gu=K-p.mat.M*(Fu+fnlj)+p.nc.sf*p.mat.Q; 