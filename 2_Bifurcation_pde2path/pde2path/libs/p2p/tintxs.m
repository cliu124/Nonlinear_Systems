function [p,t1,ts,nc]=tintxs(p,t0,ts,dt,nt,nc,pmod,smod,nffu,varargin)
% TINTXS: time integration with time-series output, LU-decomp. of M+dt*K, 
%  Here M,K pre-assembled from simple semilinear implementation.
%
%  [p,t1,ts,nc]=tintxs(p,t0,ts,dt,nt,nc,pmod,smod,nffu,varargin)
%
% t0=initial time, ts=time-series (t and residual), dt=step size, 
% nt= number of steps, nc=counter for file naming.
% plot each pmod-th step, save each smod-th step (only u, if p0 exists)
% nonlinearity f taken from nffu(p,p.u) 
% varargin: replacement stiffness matrix.
%
% Returns: t1=end time, and resulting ts, nc
%
% See also tint, tints, stanparam
r=norm(resi(p,p.u),'inf'); ts=[ts [t0;r]]; % put time and residual into ts
if ~exist([p.file.dir '/pt0.mat'],'file'); p.t=0; 
    r=norm(resi(p,p.u),'inf'); p.ts=[p.t;r]; 
    p.file.count=0; p.sol.ptype=0;  p.fuha.savefu(p); 
end
noa=nargin-9; 
if noa==0; Lam=p.mat.M+dt*p.mat.K; 
    try; if any(p.mat.Kadv) Lam=Lam+dt*p.mat.Kadv; end; end
elseif noa==1; Lam=p.mat.M+dt*varargin{1}; 
    try; if any(p.mat.Kadv) Lam=Lam+dt*p.mat.Kadv; end; end 
else  Lam=p.mat.M+dt*(varargin{1}+varargin{2});
end
%Lam=p.mat.M+dt*p.mat.K; if any(p.mat.Kadv) Lam=Lam+dt*p.mat.Kadv; end
[L,U,P,Q,R]=lu(Lam); t=t0; n=0; 
while(n<nt) % integration loop
  f=nffu(p,p.u); f=p.mat.M*f; g=p.mat.M*p.u(1:p.nu)+dt*f; 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); t=t+dt;  n=n+1;
  if(mod(n,pmod)==0); 
      r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
      tits=['t=' mat2str(t,4) ', r=' mat2str(r,3)];
      plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      title(['u_1, ' tits],'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
      drawnow; 
  end
  if(mod(n,smod)==0); 
    ps=p; p=[]; p.t=t;p.ts=ts; p.u=ps.u; 
    fname=[ps.file.pname,sprintf('%i',nc+n),'.mat']; save(fname,'p'); 
    p=ps;
  end
end 
t1=t; nc=nc+nt; 