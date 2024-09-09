function [u,res,iter,Gu,Glam,p]=nloopnleq1(p,u1)
% nloopnleq1: Newton loop by calling NLEQ1
% (done for p.sw.newt=2) 
%
[iopt,wk]=nleq1init(p); xscal=1e-3*ones(p.nu+p.nc.nq,1); rtol=p.nc.tol; 
ua=u2au(p,u1); 
% wk could be reused by storing in, e.g., p.wk, but here we refrain from this
[u2a,info,wk]=nleq1hu('nleq1f',ua,xscal,rtol,iopt,p,wk) ;
u=au2u(p,u2a); 
p.info=info; p.info=rmfield(p.info,{'xscal','xiter'}); 
r=resi(p,u); res=norm(r,p.sw.norm); % residual 
% some postprocessing
iter=info.niter; [Gu,Glam]=getder(p,u,r);
