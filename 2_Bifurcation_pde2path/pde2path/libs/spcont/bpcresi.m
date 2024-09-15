function r=bpcresi(p,u)
% bpcresi: residual for extended system for BP continuation
%
% H(U)=[G+mu*M*psi; Gu^T*psi; ||psi||^2-1; Glam^T*psi] \in R^(2*nu+2*nq+2) 
% 
% or rather 
%
% H(U)=[(G+mu*M*psi)_p; (Gu^T*psi)_p;   % PDE part
%       (G+mu*M*psi)_q; (Gu^T*psi)_q;   % q part
%       ||psi||^2-1; Glam^T*psi]        % normaliz.and constraint 
% 
% where: 
% G=regular problem part=[Gpde;q] (dimension nu+nq)
% psi=adjoint kernel eigenvector  (dimension nu+nq)
% Gu=derivative of G (incl.q) wrt (u,uaux)  (dimension nu+nq x nu+nq) 
% Glam=derivative of G wrt lam, lam='new' primary cont.par for the BPcont
% 
p=bpreduce(p); % turn off bpcont and remove new primary para.:
u1=[u(1:p.nu);u(2*p.nu+1:2*p.nu+p.naux)]; % normal pde and aux vars part
r=pderesi(p,u1); % pde part 
if(p.nc.nq>0) rq=p.fuha.qf(p,u1); else rq=[]; end;  % q part
r=[r;rq];
Glam=getGlam(p,u1,r);Gu=getGu(p,u1,r); % includes \pa_u q
if 0; jsw=p.sw.jac; p.sw.jac=0; Gu1=getGu(p,u1,r); Jd=abs(Gu1-Gu); full(Jd(1:5,1:5)), 
    full(Gu(1:5,1:5)), full(Gu1(1:5,1:5)), pause, p.sw.jac=jsw; end 
psi=[u(p.nu+1:2*p.nu); u(2*p.nu+p.naux+2:end)]; % adj.Evec
mu=u(2*p.nu+p.naux+1); % the 'zero' evalue 
M=getM(p); %M=p.mat.M(1:p.nu,1:p.nu);
M=[M                   , zeros(p.nu,p.nc.nq); 
    zeros(p.nc.nq,p.nu), zeros(p.nc.nq,p.nc.nq)];
r=r+mu*M*psi; 
rlin=Gu'*psi; 
rnorm=norm([psi; u(2*p.nu+p.naux+1:end)])^2-1; 
rcon=psi'*Glam; 
r=[r(1:p.nu); rlin(1:p.nu); rq; rlin(p.nu+1:end); rnorm; rcon];   % sort into residual 