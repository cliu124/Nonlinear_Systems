function r=fpcresi(p,u)
% fpcresi: residual for extended system for FP continuation
%
% H(U)=[G; Gu^T*phi; ||phi||^2-1]; \in R^(2*nu+2*nq+1) 
% 
% or rather 
%
% H(U)=[(G)_p; (Gu*phi)_p; % PDE part
%        q; (Gu^T*phi)_q;   % q part
%       ||phi||^2-1]        % normaliz. 
% 
% where: 
% G=regular problem part=[Gpde;q] (dimension nu+nq)
% phi=kernel eigenvector  (dimension nu+nq)
% Gu=derivative of G (incl.q) wrt (u,uaux)  (dimension nu+nq x nu+nq) 
% 
p=spreduce(p); % turn off spcont and remove new foco primary param 
u1=[u(1:p.nu); u(2*p.nu+1:2*p.nu+p.naux)]; % normal pde part and pars 
r=pderesi(p,u1); if(p.nc.nq>0) rq=p.fuha.qf(p,u1); else rq=[]; end % pde and aux part:
% add linearization residual: Gu*phi, phi=norm. eigenvector in kernel 
Gu=getGu(p,u1,[r;rq]); % includes \pa_u q
phi=[u(p.nu+1:2*p.nu);u(2*p.nu+p.naux+1:end)]; rlin=Gu*phi;
r=[r; rlin(1:p.nu); rq; rlin(p.nu+1:end)]; % sort into residual 
r=[r; norm(phi)^2-1];  % append normalization of eigenfunction <phi,phi>=1:  