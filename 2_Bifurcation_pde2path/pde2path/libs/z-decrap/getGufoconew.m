function Gu=getGufoco(p,u,r)
% getGufoco: FP continuation extension of Gu. The assumption is that initially 
% this will be called after spcontini, so that the old prim. param. is at p.nc.ilam(2).
%
%  v=upde normal unknowns part,   w=original active uaux part
%  ph=u(p.nu/2+1:p.nu) pde part of kernel eigenvector
%  om=u(p.nu+p.naux+1:end) active non-primary uaux part of kernel eigenvector 
%  G=pde, Q=auxiliary equations
%  Gv=partial derivative G_v, others analogous,  Gvvph=\pa_v(G_v*ph)
p=spreduce(p); % set regular case sizes to compute G and phi
u1=[u(1:p.nu); u(2*p.nu+1:2*p.nu+p.naux)]; % org variables (no phi)
nu1=length(u1); % this is p.nu+p.naux (with p.nu halved)
r1=[r(1:p.nu); r(2*p.nu+1:2*p.nu+p.nc.nq)]; % org variables
ph=u(p.nu+1:2*p.nu); om=u(2*p.nu+p.naux+1:end); % eigenvector (pde and aux part) 
Gv=getGupde(p,u1,r1);     % get analytic pde part linearization
Gvph=Gv*u(p.nu+1:2*p.nu); % take reference from residual
if p.sw.spjac==1; Gvvph=p.fuha.spjac(p,u); % user def. 
else; Gvvph=sparse(p.nu,p.nu); % home-made expensive way 
  for j=1:p.nu  
   Gv1=getGupde(p,u1+p.nc.del*ej(j,nu1),r1); % perturbed pde-part lin.
   r2=(Gv1*ph-Gvph)/p.nc.del;
   % add column with finite diff. approx.:
   Gvvph=Gvvph+sparse(1:p.nu,j,r2,p.nu,p.nu); % Gvvph(:,j-1)=(r1-r)/del; 
  end   
end
% finite difference for mixed pde-aux part \pa_w(G_v*ph) 
Gvwph=sparse(p.nu,p.nc.nq+1); 
% first perturb by old primary parameter (other colums (nq>0) later) 
Gv1=getGupde(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(1),nu1),r1); 
% part of mixed deriv. multiplied with ph (pde part of evec.)
r1=(Gv1*ph-Gvph)/p.nc.del; 
Gvwph=Gvwph+sparse(1:p.nu,1,r1,p.nu,p.nc.nq+1);
if(p.nc.nq==0)
  r2=pderesi(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(1),nu1));
  Gw=(r2-r(1:p.nu))/p.nc.del;  % Gw finite diff via pderesi:
  % assemble Gu from block matrices; last row is the normalization
  Gu=[[ Gv   0*Gv   Gw  ];  % pa_(u,phi,w) G
      [Gvvph   Gv   Gvwph]; % pa_(u,phi,w) (G_u*phi)
      [0*ph' p.phi0'   0  ]];
  return; 
end
% now case nq>0; exclude new primary parameter (relevant for G_lam), and initially also old one:
Gwvom=sparse(p.nu,p.nu);   
for k=1:p.nc.nq 
       Gv1=getGupde(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(k+1),nu1),r1); 
       % part of mixed deriv. multiplied with phi (pde part of evec.)
       r2=(Gv1*ph-Gvph)/p.nc.del;
       Gvwph=Gvwph+sparse(1:p.nu,k+1,r2,p.nu,p.nc.nq+1);
       % part of mixed deriv. multiplied with components of om (aux part of evec.)
       Gwvom=Gwvom+u(2*p.nu+p.naux+k)*(Gv1-Gv)/p.nc.del;
end
% derivatives of q_v*om: 
Qvvph=sparse(p.nc.nq,p.nu);      % \pa_v(q_v*phi)
Qwvom=sparse(p.nc.nq,p.nu);      % \pa_u(q_w*om) 
Qvwph=sparse(p.nc.nq,p.nc.nq+1); % \pa_w(q_v*phi) 
rq=r(2*p.nu+1:2*p.nu+p.nc.nq);   % aux. eqn. residual
%Qvph=p.fuha.qfder(p,u1)*ph;     % q_v*phi, aux. eqn. linearization residual
Qvph=r(2*p.nu+p.nc.nq+1:length(r)-1)'; % same, but a poor iguess might cause divergence
if p.sw.qjac==1 % analytical q_u 
   Qv=p.fuha.qfder(p,u1); spqjac=0; if isfield(p.sw,'spqjac'); spqjac=p.sw.spqjac; end 
   if spqjac==1; Qvvph=p.fuha.spqjac(p,u); % analytical pa_u(q_u*phi) 
   else
     for j=1:p.nu
       Qv1=p.fuha.qfder(p,u1+p.nc.del*ej(j,nu1)); % perturbed pde-part lin.
       r2=(Qv1*ph-Qvph)/p.nc.del;
       Qvvph=Qvvph+sparse(1:p.nc.nq,j,r2,p.nc.nq,p.nu);
     end
   end
   Qv1=p.fuha.qfder(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(1),nu1));   % \pa_w(q_u*phi) 
   r2=(Qv1*ph-Qvph)/p.nc.del;
   Qvwph=Qvwph+sparse(1:p.nc.nq,1,r2,p.nc.nq,p.nc.nq+1);
   for k=1:p.nc.nq
       Qv1=p.fuha.qfder(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(k+1),nu1)); 
       r2=(Qv1*ph-Qvph)/p.nc.del;
       Qvwph=Qvwph+sparse(1:p.nc.nq,k+1,r2,p.nc.nq,p.nc.nq+1);
       Qwvom=Qwvom+u(2*p.nu+p.naux+k)*(Qv1-Qv)/p.nc.del;
   end
else % all second order numerically
   Qv=sparse(p.nc.nq,p.nu); 
   for j=1:p.nu
       r1=p.fuha.qf(p,u1+p.nc.del*ej(j,nu1));
       Qv=Qv+sparse(1:p.nc.nq,j,(r1-rq)/p.nc.del,p.nc.nq,p.nu);           
       for jj=1:p.nu  % second derivatives: SLOW!
          rpp=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)+ej(jj,nu1))); 
          rmm=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)+ej(jj,nu1))); 
          rpm=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)-ej(jj,nu1)));
          rmp=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)-ej(jj,nu1)));
          r2=u(p.nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qvvph=Qvvph+sparse(1:p.nc.nq,j,r2,p.nc.nq,p.nu);
       end
   %   mixed derivatives, first (initial run) wrt the old primary parameter
       kk=p.nu+p.nc.ilam(1);
       rpp=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
       rmm=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
       rpm=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
       rmp=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
       r2=u(p.nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
       Qvwph=Qvwph+sparse(1:p.nc.nq,1,r2,p.nc.nq,p.nc.nq+1);
       %   now with respect to other old active aux. vars
       for k=2:p.nc.nq+1 
          kk=p.nu+p.nc.ilam(k);
          rpp=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
          rmm=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
          rpm=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
          rmp=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
          r2=u(2*p.nu+p.naux+k-1)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qwvom=Qwvom+sparse(1:p.nc.nq,j,r2,p.nc.nq,p.nu);
          r2=u(p.nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qvwph=Qvwph+sparse(1:p.nc.nq,k-1,r2,p.nc.nq,p.nc.nq+1);
       end
   end
end % p.sw.qjac=1
% Now build G_U (U=(u,phi,w), i.e. everything except prim parameter
% Notation: 'lam' here means the old primary parameter
% In case p.nc.nq=0 above this was named 'Gw' while here Gw stands for 
% other aux. vars. derivatives.
Gwlam=sparse(p.nu,p.nc.nq+1); %pa_(lam,w) G
Gwwom=sparse(p.nu,p.nc.nq+1); %pa_(lam,om)(G_u*[phi;om])
Qwlam=sparse(p.nc.nq,p.nc.nq+1); Qwwom=sparse(p.nc.nq,p.nc.nq+1);
for k=2:p.nc.nq+2
   k1=p.nu+p.nc.ilam(k-1);
   % Gw finite diff via pde:
   r1=pderesi(p,u1+p.nc.del*ej(k1,nu1));
   Gwlam=Gwlam+sparse(1:p.nu,k-1,(r1-r(1:p.nu))/p.nc.del,p.nu,p.nc.nq+1);           
   % Qw finite diff:
   r1=p.fuha.qf(p,u1+p.nc.del*ej(k1,nu1));
   Qwlam=Qwlam+sparse(1:p.nc.nq,k-1,(r1-rq)/p.nc.del,p.nc.nq,p.nc.nq+1);
   % second derivatives:
   for kk=3:p.nc.nq+2 % start at 2 to skip primary parameter
      kk1=p.nu+p.nc.ilam(kk-1);
      rpp=p.fuha.qf(p,u1+p.nc.del*(ej(k1,nu1)+ej(kk1,nu1))); rmm=p.fuha.qf(p,u1-p.nc.del*(ej(k1,nu1)+ej(kk1,nu1)));
      rpm=p.fuha.qf(p,u1+p.nc.del*(ej(k1,nu1)-ej(kk1,nu1))); rmp=p.fuha.qf(p,u1-p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)));
      r2=u(2*p.nu+p.naux+kk-2)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Qwwom=Qwwom+sparse(1:p.nc.nq,k-1,r2,p.nc.nq,p.nc.nq+1);
      % now for Gww... 
      rpp=pderesi(p,u1+p.nc.del*(ej(k1,nu1)+ej(kk1,nu1))); rmm=pderesi(p,u1-p.nc.del*(ej(k1,nu1)+ej(kk1,nu1))); 
      rpm=pderesi(p,u1+p.nc.del*(ej(k1,nu1)-ej(kk1,nu1))); rmp=pderesi(p,u1-p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)));
      r2=u(2*p.nu+p.naux+kk-2)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Gwwom=Gwwom+sparse(1:p.nu,k-1,r2,p.nu,p.nc.nq+1);
   end
end
% remove spcont-parameter (added later via Jac=[Gu Glam; tau]) 
Gw=Gwlam(:,2:end); Qw=Qwlam(:,2:end); 
% assemble Gu from block matrices; last row is the normalization
zeroGv=sparse(p.nu,p.nu); zeroGw=sparse(p.nu,p.nc.nq); 
zeroQw=sparse(p.nc.nq,p.nc.nq); zeroQv=sparse(p.nc.nq,p.nu);
zerovt=sparse(1,p.nu);    zerowt=sparse(1,p.nc.nq+1);
% Gv=nu x nu, Gwlam=nu x nq+1, Gw=nu x nq
Gu=[[    Gv      zeroGv    Gwlam     zeroGw];  % (pa_v, pa_phi, pa_w, pa_om) G 
    [Gvvph+Gwvom   Gv   Gvwph+Gwwom    Gw  ];  % (pa_v .. pa_om) (G_U*(phi,om))                         
    [    Qv      zeroQv    Qwlam     zeroQw];  % (pa_v .. pa_om) Q
    [Qvvph+Qwvom   Qv   Qvwph+Qwwom    Qw  ];  % (pa_v .. pa_om) (Q_U*(phi,om))
    [  zerovt     2*ph'    zerowt     2*om']];
end 