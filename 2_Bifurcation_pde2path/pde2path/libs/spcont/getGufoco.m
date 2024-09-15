function Gu=getGufoco(p,u,r)
% getGufoco: FP continuation Jacobian 
% see fpcresi for the pertinent extended system H 
%     spcontini for the setup of variables; 
% roughly: old prim. param. is at p.nc.ilam(2).
%
%  v=upde normal unknowns part,   w=original active uaux part
%  ph=u(p.nu/2+1:p.nu) pde part of kernel eigenvector
%  om=u(p.nu+p.naux+1:end) active non-primary uaux part of kernel eigenvector 
%  G=pde, Q=auxiliary equations
%  Gv=partial derivative G_v, others analogous,  Gvvph=\pa_v(G_v*ph)
jtest=0; 
p=spreduce(p); % set regular case sizes to compute G and phi
nu=p.nu; nq=p.nc.nq; 
u1=[u(1:nu); u(2*nu+1:2*nu+p.naux)]; % org variables (no phi)
nu1=length(u1); % this is nu+p.naux (with nu halved)
r1=[r(1:nu); r(2*nu+1:2*nu+nq)]; % org variables
%r1b=resi(p,u1); max(abs(r1-r1b)), pause 
ph=u(nu+1:2*nu); om=u(2*nu+p.naux+1:end); % eigenvector (pde and aux part) 
%jaccheck(p,0.1); pause 
Gv=getGupde(p,u1,r1); % pde-part linearization
Gvph=Gv*u(nu+1:2*nu); % pde-part of Gu*phi
switch p.sw.spjac
    case 1; Gvvph=p.fuha.spjac(p,u); % user def. 
    case 0; Gvvph=getGuphiu(p,u1,ph); 
    otherwise; Gvvph=sparse(nu,nu); % home-made expensive way 
  for j=1:nu  % loop over v 
      if 1 % forward 
   up=u1+p.nc.del*ej(j,nu1); r2=pderesi(p,up); 
   Gv1=getGupde(p,up,r2); % perturbed pde-part lin.
   r3=(Gv1*ph-Gvph)/p.nc.del;   
      else % backward 
   up=u1-p.nc.del*ej(j,nu1); r2=pderesi(p,up); 
   Gv1=getGupde(p,up,r2); % perturbed pde-part lin.
   r3=-(Gv1*ph-Gvph)/p.nc.del;        
      end
   Gvvph=Gvvph+sparse(1:nu,j,r3,nu,nu);    
  end   
end
if jtest; Gvvph1=getGuphiu(p,u1,ph); Jd=abs(Gvvph1-Gvvph); 
    full(Gvvph(1:5,1:5)), full(Gvvph1(1:5,1:5)),max(max(Jd)), pause; end 
% finite difference for mixed pde-aux part \pa_w(G_v*ph) 
Gvwph=sparse(nu,nq+1); 
% first perturb by old primary parameter (other colums (nq>0) later) 
up=u1+p.nc.del*ej(nu+p.nc.ilam(1),nu1); r1=resi(p,up); %up(1:5)', r1(1:5)', pause 
Gv1=getGupde(p,up,r1); 
if 0
jsw=p.sw.jac; p.sw.jac=1; Gv1b=getGupde(p,up,r1);  Gvd=abs(Gv1b-Gv1); 
full(Gv1(1:5,1:5)), full(Gv1b(1:5,1:5)), full(Gvd(1:5,1:5)), pause; p.sw.jac=jsw; 
end 
% part of mixed deriv. multiplied with ph (pde part of evec.)
r1=(Gv1*ph-Gvph)/p.nc.del; 
Gvwph=Gvwph+sparse(1:nu,1,r1,nu,nq+1); % 1 column, der. wrt old primary 
if(nq==0)
  r2=pderesi(p,u1+p.nc.del*ej(nu+p.nc.ilam(1),nu1));
  Gw=(r2-r(1:nu))/p.nc.del;  % Gw finite diff via pderesi:
  % assemble Gu from block matrices; last row is the normalization
  Gu=[[ Gv   0*Gv   Gw  ];  % pa_(u,phi,w) G
      [Gvvph   Gv   Gvwph]; % pa_(u,phi,w) (G_u*phi)
      [0*ph' 2*ph'   0  ]];
  return; 
end
% now case nq>0; exclude new primary parameter (relevant for G_lam), and 
% initially also old one:
try quuphisw=p.sw.quuphi; catch quuphisw=1; end 
Gwvom=sparse(nu,nu);   
for k=1:nq 
       up=u1+p.nc.del*ej(nu+p.nc.ilam(k+1),nu1); r2=pderesi(p,up); 
       Gv1=getGupde(p,up,r2); 
       % part of mixed deriv. multiplied with phi (pde part of evec.)
       r3=(Gv1*ph-Gvph)/p.nc.del;
       Gvwph=Gvwph+sparse(1:nu,k+1,r3,nu,nq+1);
       % part of mixed deriv. multiplied with components of om (aux part of evec.)
       Gwvom=Gwvom+u(2*nu+p.naux+k)*(Gv1-Gv)/p.nc.del;
end
% derivatives of q_v*om: 
Qwvom=sparse(nq,nu);   % \pa_u(q_w*om) 
Qvwph=sparse(nq,nq+1); % \pa_w(q_v*phi) 
rq=r(2*nu+1:2*nu+nq);  % q-residual
Qv=p.fuha.qfder(p,u1); % pa_u q, see below if this is not provided! 
Qvph=Qv*ph;      % q_v*phi, aux. eqn. linearization residual
spqjac=1; if isfield(p.sw,'spqjac'); spqjac=p.sw.spqjac; end 
if spqjac==1; 
    Qvvph=quuphi(p,u); % analytical pa_u(q_u*phi) (usually zero)    
else; Qvvph=sparse(nq,nu); % \pa_v(q_v*phi) % numerical,very expensive! 
  for j=1:nu
     up=u1+p.nc.del*ej(nu+p.nc.ilam(k+1),nu1);
     Qv1=p.fuha.qfder(p,up); % perturbed pde-part lin.
     r2=(Qv1*ph-Qvph)/p.nc.del;
     Qvvph=Qvvph+sparse(1:nq,j,r2,nq,nu);
  end
end % Qvvph
up=u1+p.nc.del*ej(nu+p.nc.ilam(1),nu1); 
Qv1=p.fuha.qfder(p,up);   % \pa_w(q_u*phi) 
r2=(Qv1*ph-Qvph)/p.nc.del; 
Qvwph=Qvwph+sparse(1:nq,1:nq,r2,nq,nq+1);
for k=1:nq
    up=u1+p.nc.del*ej(nu+p.nc.ilam(k+1),nu1); 
    Qv1=p.fuha.qfder(p,up); 
    r2=(Qv1*ph-Qvph)/p.nc.del;
    Qvwph=Qvwph+sparse(1:nq,k+1,r2,nq,nq+1);
    Qwvom=Qwvom+u(2*nu+p.naux+k)*(Qv1-Qv)/p.nc.del;
end
% Now build G_U (U=(u,phi,w), i.e. everything except prim parameter
% Notation: 'lam' here means the old primary parameter
% In case nq=0 above this was named 'Gw' while here Gw stands for 
% other aux. vars. derivatives.
Gwlam=sparse(nu,nq+1); %pa_(lam,w) G
Gwwom=sparse(nu,nq+1); %pa_(lam,om)(G_u*[phi;om])
Qwlam=sparse(nq,nq+1); %pa_(lam,w) Q
Qwwom=sparse(nq,nq+1); %pa_(w,om) Q 
for k=2:nq+2
   k1=nu+p.nc.ilam(k-1);
   % Gw finite diff via pde:
   up=u1+p.nc.del*ej(k1,nu1); 
   r1=pderesi(p,up); 
   Gwlam=Gwlam+sparse(1:nu,k-1,(r1-r(1:nu))/p.nc.del,nu,nq+1);           
   % Qw finite diff:
   r1=p.fuha.qf(p,up); 
   Qwlam=Qwlam+sparse(1:nq,k-1,(r1-rq)/p.nc.del,nq,nq+1);
   % second derivatives:
   for kk=3:nq+2 % start at 2 to skip primary parameter
      kk1=nu+p.nc.ilam(kk-1);
      rpp=p.fuha.qf(p,u1+p.nc.del*(ej(k1,nu1)+ej(kk1,nu1))); 
      rmm=p.fuha.qf(p,u1-p.nc.del*(ej(k1,nu1)+ej(kk1,nu1)));
      rpm=p.fuha.qf(p,u1+p.nc.del*(ej(k1,nu1)-ej(kk1,nu1))); 
      rmp=p.fuha.qf(p,u1-p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)));
      r2=u(2*nu+p.naux+kk-2)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Qwwom=Qwwom+sparse(1:nq,k-1,r2,nq,nq+1);
      % now for Gww... 
      rpp=pderesi(p,u1+p.nc.del*(ej(k1,nu1)+ej(kk1,nu1))); 
      rmm=pderesi(p,u1-p.nc.del*(ej(k1,nu1)+ej(kk1,nu1))); 
      rpm=pderesi(p,u1+p.nc.del*(ej(k1,nu1)-ej(kk1,nu1))); 
      rmp=pderesi(p,u1-p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)));
      r2=u(2*nu+p.naux+kk-2)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Gwwom=Gwwom+sparse(1:nu,k-1,r2,nu,nq+1);
   end
end
% remove spcont-parameter (added later via Jac=[Gu Glam; tau]) 
Gw=Gwlam(:,2:end); Qw=Qwlam(:,2:end); 
% assemble Gu from block matrices; last row is the normalization
zeroGv=sparse(nu,nu); zeroGw=sparse(nu,nq); 
zeroQw=sparse(nq,nq); zeroQv=sparse(nq,nu);
zerovt=sparse(1,nu);    zerowt=sparse(1,nq+1);
% Gv=nu x nu, Gwlam=nu x nq+1, Gw=nu x nq
%        v      phi-pde      w     om=phi-q
Gu=[[    Gv      zeroGv    Gwlam     zeroGw];  % (pa_v, pa_phi, pa_w, pa_om) G 
    [Gvvph+Gwvom   Gv   Gvwph+Gwwom    Gw  ];  % (pa_v .. pa_om) (G_U*(phi,om))                         
    [    Qv      zeroQv    Qwlam     zeroQw];  % (pa_v .. pa_om) Q
    [Qvvph+Qwvom   Qv   Qvwph+Qwwom    Qw  ];  % (pa_v .. pa_om) (Q_U*(phi,om))
    [  zerovt     2*ph'    zerowt     2*om']]; % der. of normalization ||phi||^2-1
return
% -------------------------------------------------------------------------
% rem: use this if no qfder is provided (extremely inefficient!) 
rq=r(2*nu+1:2*nu+nq);   Qv=sparse(nq,nu); 
for j=1:nu
   r1=p.fuha.qf(p,u1+p.nc.del*ej(j,nu1));
   Qv=Qv+sparse(1:nq,j,(r1-rq)/p.nc.del,nq,nu);           
   for jj=1:nu  % second derivatives: SLOW!
      rpp=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)+ej(jj,nu1))); 
      rmm=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)+ej(jj,nu1))); 
      rpm=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)-ej(jj,nu1)));
      rmp=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)-ej(jj,nu1)));
      r2=u(nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Qvvph=Qvvph+sparse(1:nq,j,r2,nq,nu);
   end
%   mixed derivatives, first (initial run) wrt the old primary parameter
   kk=nu+p.nc.ilam(1);
   rpp=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
   rmm=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
   rpm=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
   rmp=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
   r2=u(nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
   Qvwph=Qvwph+sparse(1:nq,1,r2,nq,nq+1);
   %   now with respect to other old active aux. vars
   for k=2:nq+1 
      kk=nu+p.nc.ilam(k);
      rpp=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
      rmm=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)+ej(kk,nu1))); 
      rpm=p.fuha.qf(p,u1+p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
      rmp=p.fuha.qf(p,u1-p.nc.del*(ej(j,nu1)-ej(kk,nu1)));
      r2=u(2*nu+p.naux+k-1)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Qwvom=Qwvom+sparse(1:nq,j,r2,nq,nu);
      r2=u(nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
      Qvwph=Qvwph+sparse(1:nq,k-1,r2,nq,nq+1);
   end
end