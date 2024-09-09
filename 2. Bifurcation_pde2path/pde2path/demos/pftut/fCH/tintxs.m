function [p,t1,ts,nc]=tintxs(p,t0,ts,dt,nt,nc,pmod,smod,nffu,varargin)
% Mod of tintxs for fCH since stiffness matrix assembled in non-standard way
r=norm(resi(p,p.u),'inf'); ts=[ts [t0;r]]; % put time and residual into ts
if ~exist([p.file.dir '/pt0.mat'],'file'); p.t=0; 
    r=norm(resi(p,p.u),'inf'); p.ts=[p.t;r]; 
    p.file.count=0; p.sol.ptype=0;  p.fuha.savefu(p); 
end
par=p.u(p.nu+1:end); epsi=par(3); eps2=epsi^2; 
Ms=p.mat.Ms; Mnl=[[Ms 0*Ms];[0*Ms Ms]];
Ks=p.mat.Ks; Ksys=[[0*Ks -eps2*Ks]; [eps2*Ks Ms]]; 
Lam=p.mat.M+dt*Ksys; % spy(Ksys); pause; spy(Mnl); 
[L,U,P,Q,R]=lu(Lam); t=t0; n=0; 
while(n<nt) % integration loop
  f=nffu(p,p.u); f=Mnl*f; g=p.mat.M*p.u(1:p.nu)+dt*f; 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); t=t+dt;  n=n+1;
  if(mod(n,pmod)==0); 
      r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
      tits=['t=' mat2str(t,4) ', r=' mat2str(r,3)];
      plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      %title(['u_1, ' tits],'fontsize',p.plot.fs); 
      title(['t=' mat2str(t,4) ], 'fontsize',14); 
      set(gca,'FontSize',p.plot.fs); 
      drawnow; 
  end
  if(mod(n,smod)==0); 
    ps=p; p=[]; p.t=t;p.ts=ts; p.u=ps.u; 
    fname=[ps.file.pname,sprintf('%i',nc+n),'.mat']; save(fname,'p'); 
    p=ps;
  end
end 
t1=t; nc=nc+nt; 