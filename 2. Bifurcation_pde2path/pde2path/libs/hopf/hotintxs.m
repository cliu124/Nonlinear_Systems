function [p,t1,ts,nc]=hotintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,nffu,pti)
% HOTINTXS: time integration with LU-decomp. of M+dt*K, 
%  adapted to Hopf; ts-output of || u0-u(nT) || 
%  Here M,K pre-assembled from simple semilinear implementation.
%
%  [p,t1,ts,nc]=tintxs(p,u0,t0,ts,npp,nt,nc,tsmod,pmod,smod,nffu,pti) 
%
% t0=initial time, ts=time-series, npp=# points per period 
% nt= number of steps, nc=counter for u output, 
% plot each pmod-th step, save each smod-th step (only u, if p0 exists)
% ts-output each tsmod step
% pti=point indices for values to be put in ts
% 
% nonlinearity f taken from nffu(p,p.u) 
%
% Returns: t1=end time, and resulting ts, nc
%
% See also tint, tints, stanparam
vals=p.u(pti); nd=norm(p.u(1:p.nu)-u0(1:p.nu),'inf'); 
ts=[ts [t0; vals; nd]]; % put time and values into ts
if ~exist([p.file.dir '/pt0.mat'],'file'); p.t=0; 
    p.file.count=0; p.sol.ptype=0;  p.fuha.savefu(p); 
end
dt=p.hopf.T/npp; % npp=# points per period 
Lam=p.mat.M+dt*p.mat.K;
[L,U,P,Q,R]=lu(Lam); t=t0; n=0; 
while(n<nt) % integration loop
  f=nffu(p,p.u); f=p.mat.M*f; g=p.mat.M*p.u(1:p.nu)+dt*f; 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); t=t+dt;  n=n+1;
  if(mod(n,tsmod)==0); % put time and values into ts
       vals=p.u(pti); nd=norm(p.u(1:p.nu)-u0(1:p.nu),'inf'); 
       ts=[ts [t; vals; nd]]; end 
  if(mod(n,pmod)==0); % plot 
      tits=['t=' mat2str(t,4)];
      plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      title(['u_1, ' tits],'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
      drawnow;   end
  if(mod(n,smod)==0); % save 
    ps=p; p=[]; p.t=t;p.ts=ts; p.u=ps.u; 
    fname=[ps.file.pname,sprintf('%i',nc+n),'.mat']; save(fname,'p'); 
    p=ps;
  end
end 
t1=t; nc=nc+nt; 