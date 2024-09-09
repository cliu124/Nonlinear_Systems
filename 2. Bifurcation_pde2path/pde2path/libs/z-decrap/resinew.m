function r=resi(p,u)
% RESI: return residual for PDE and auxiliary equations
%
%  r=resi(p,u)
%
% See also pderesi, getGu
if (p.sw.spcont==0)          % regular continuation 
    r=pderesi(p,u);           % pde part
  if(p.nc.nq>0) r=[r;p.fuha.qf(p,u)]; end  % possibly add auxiliary functions
end
if p.sw.spcont==1 % BP continuation 
    ilam=p.nc.ilam; p=bpreduce(p); % turn off spectral cont and remove new fold/branch pt primary para.:
    u1=[u(1:p.nu);u(2*p.nu+1:2*p.nu+p.naux)]; % normal pde part
    r=pderesi(p,u1); % pde part 
    if(p.nc.nq>0) rq=p.fuha.qf(p,u1); else rq=[]; end
    Glam=getGlam(p,u1,r);  Gu=getGu(p,u1,[r;rq]); 
    psi=u(p.nu+1:2*p.nu); mu=u(2*p.nu+ilam(end)); % the (zero) eigenvalue/vector 
    r=r+mu*p.mat.M(1:p.nu,1:p.nu)*psi; rlin=Gu'*psi; 
    rnorm=norm([psi; 0*u(2*p.nu+p.naux+1:end)])^2-1; rcon=psi'*Glam; 
    r=[r; rlin(1:p.nu); rq; rlin(p.nu+1:end); rnorm; rcon];   % sort into residual 
end
if p.sw.spcont==2 % fold point continuation 
    p=spreduce(p); % turn off spcont and remove new foco primary param 
    u1=[u(1:p.nu); u(2*p.nu+1:2*p.nu+p.naux)]; % normal pde part and pars 
    r=pderesi(p,u1); if(p.nc.nq>0) rq=p.fuha.qf(p,u1); else rq=[]; end % pde and aux part:
    % add linearization residual: Gu*phi, phi=norm. eigenvector in kernel 
    Gu=getGu(p,u1,[r;rq]); phi=[u(p.nu+1:2*p.nu);u(2*p.nu+p.naux+1:end)]; rlin=Gu*phi;
    r=[r; rlin(1:p.nu); rq; rlin(p.nu+1:end)]; % sort into residual 
    r=[r; p.phi0'*phi-1];  % append normalization of eigenfunction <phi,phi>=1:
  %  norm(r,'inf'),norm(rlin,'inf'),norm(p.phi0'*phi-1,'inf'), pause
end
if p.sw.spcont==3 % HP continuation 
    p=hpreduce(p); % turn off spcont and remove new foco primary param 
    u1=[u(1:p.nu); u(3*p.nu+2:3*p.nu+1+p.naux)]; % normal pde part and pars 
    r1=pderesi(p,u1); if(p.nc.nq>0) rq=p.fuha.qf(p,u1); else rq=[]; end % pde and aux part:
    % append further resis: Gu*phi_r+om*M*phr; Gu*phi-om*M*phr,
    %  phr^2+phi^2=1; 2*phr*phi=0 
    Gu=getGu(p,u1,[r1;rq]); M=p.mat.M(1:p.nu,1:p.nu); 
    phr=u(p.nu+1:2*p.nu); phi=u(2*p.nu+1:3*p.nu); om=u(3*p.nu+1); 
    r2=Gu*phr+om*M*phi; r3=Gu*phi-om*M*phr; r4=p.c*phr-1; r5=p.c*phi; 
   % fprintf('r1=%g, r2=%g, r3=%g, r4=%g, r5=%g \n',norm(r1,'inf'),norm(r2,'inf'),norm(r3,'inf'),norm(r4,'inf'),norm(r5,'inf')); 
    r=[r1; r2; r3; r4; r5]; % sort into residual
end