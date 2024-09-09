function p=tint(p,dt,nt,pmod)
% TINT: time integration, semi-implicit (Euler)steps, full FEM assembling
%
%  p=tint(p,dt,nt,pmod)
%
% See also tints, tintx, tintxs.
p.sw.sfem=0; fprintf('running tint, see figure %g\n',p.plot.ifig); 
start=tic;
a=1;c=0;f=zeros(p.nc.neq,1); % assemble mass-matrix M 
bc=p.fuha.bc(p,p.u); upde=p.mat.fill*p.u(1:p.nu); 
[K,M,F,Q,G,H,R]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,c,a,f,upde); 
M=filltrafo(p,M); % adapt in case of periodic domain
nc=0; 
while(nc<nt) 
  [c,a,f,b]=p.fuha.G(p,p.u); bc=p.fuha.bc(p,p.u); 
  upde=p.mat.fill*p.u(1:p.nu); 
  [K,F]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,c,a,f,upde);
  if(any(b)) Kadv=assemadv(p.mesh.p,p.mesh.t,b); K=K-Kadv; end 
  F=p.mat.fill'*F; K=filltrafo(p,K); % adapt in case of periodic domain
  L=M+dt*K; p.u(1:p.nu)=L\(M*p.u(1:p.nu)+dt*F); nc=nc+1;
  if(mod(nc,pmod)==0); plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      r=norm(resi(p,p.u),p.sw.norm); fprintf('time=%g, res=%g\n',dt*nc,r); 
      drawnow; end 
end 
fprintf('Timing: %g\n',toc(start));