%%
% TINTXS: time integration with time-series output, LU-decomp. of M+dt*K, 
%  Here M,K pre-assembled from simple semilinear implementation.
%
%  mod of tintxs for KSpbc
function [p,t1,ts,nc]=tintxs(p,t0,ts,dt,nt,nc,pmod,smod)
r=norm(resi(p,p.u),'inf'); ts=[ts [t0;r]]; % put time and residual into ts
Ks=p.mat.K; M=p.mat.Ms; Kx=p.mat.Kx; 
par=p.u(p.nu+1:end); al=par(1); eps=par(3); s=par(4); 
n=p.nu/2; K=[[-Ks-s*Kx -al*Ks];[-Ks -M]]; 
Lam=p.mat.M+dt*K; 
[L,U,P,Q,R]=lu(Lam); t=t0; m=0; 
while(m<nt) % integration loop
  u1=p.u(1:n); F=[0.5*Kx*(u1.^2)+eps; zeros(n,1)]; 
  g=p.mat.M*p.u(1:p.nu)+dt*F; 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); t=t+dt;  m=m+1;
  if(mod(m,pmod)==0); 
      r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
      tits=['t=' mat2str(t,4) ', r=' mat2str(r,3)];
      plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      title(['u_1, ' tits],'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
      drawnow;
  end
  if 0
  if(mod(m,smod)==0); 
    ps=p; p=[]; p.t=t;p.ts=ts; p.u=ps.u; 
    fname=[ps.file.pname,sprintf('%i',nc+n),'.mat']; save(fname,'p'); 
    p=ps;
  end
  end
end 
t1=t; nc=nc+nt; 