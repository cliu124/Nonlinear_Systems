function p=tints(p,dt,nt,pmod,nffu,varargin)
% TINTS: time integration with LU-decomp. of M+dt*K, M,K pre-assembled 
%  do nt steps; from simple semilinear implementation.
%
%   [p]=tints(p,dt,nt,pmod,nffu,varargin)
%
% * dt=step size, nt= number of steps, plot each pmod-th step, 
%   nonlinearity f taken from nffu(p,p.u) 
% * varargin: replacement stiffness matrix.
%
% See also tint, tintxs, stanparam
noa=nargin-5; 
if noa==0; Lam=p.mat.M+dt*p.mat.K; if any(p.mat.Kadv) Lam=Lam+dt*p.mat.Kadv; end
elseif noa==1; Lam=p.mat.M+dt*varargin{1}; if any(p.mat.Kadv) Lam=Lam+dt*p.mat.Kadv; end
else  Lam=p.mat.M+dt*(varargin{1}+varargin{2});
end
[L,U,P,Q,R]=lu(Lam); n=0; 
while(n<nt) % integration loop
  f=nffu(p,p.u); g=p.mat.M*p.u(1:p.nu)+dt*p.mat.M*f; 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); n=n+1; 
  if(mod(n,pmod)==0); 
   plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
  end
end 