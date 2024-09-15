function [p,t1,ts,nc]=tintx(p,t0,ts,dt,nt,nc,pmod,smod)
% TINTX: time integration with time series output, full G assembly
%
%  [p,t1,ts,nc] = tintx(p,t0,ts,dt,nt,nc,pmod,smod)
%
% See also tintxs, tint, tints
p.sw.sfem=0; t=t0; r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
if ~exist([p.file.dir '/pt0.mat'],'file'); p.t=0; 
    r=norm(resi(p,p.u),'inf'); p.ts=[p.t;r]; 
    p.file.count=0; p.sol.ptype=0;  p.fuha.savefu(p); 
end
a=1;c=0;f=zeros(p.nc.neq,1); bc=p.fuha.bc(p,p.u); % assemble mass-matrix M 
upde=p.mat.fill*p.u(1:p.nu); 
[K,M,F,Q,G,H,R]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,c,a,f,upde); 
M=filltrafo(p,M); % adapt in case of periodic domain
n=0; 
while(n<nt) 
  [c,a,f,b]=p.fuha.G(p,p.u); bc=p.fuha.bc(p,p.u); upde=p.mat.fill*p.u(1:p.nu); 
  [K,F]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,c,a,f,upde);
  if(any(b)) Kadv=assemadv(p.mesh.p,p.mesh.t,b); K=K-Kadv; end 
  F=p.mat.fill'*F; K=filltrafo(p,K); % adapt in case of periodic domain
  L=M+dt*K; p.u(1:p.nu)=L\(M*p.u(1:p.nu)+dt*F); t=t+dt; n=n+1;
  r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
  if(mod(n,pmod)==0); 
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