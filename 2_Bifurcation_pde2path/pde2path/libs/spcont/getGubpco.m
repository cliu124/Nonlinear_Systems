function Gu=getGubpco(p,u,r)
% getGubpco: BP continuation jacobian  
%
% see bpcresi for the pertinent extended system H 
%     bpcontini for the setup of variables; 
% roughly: old prim. param. is at p.nc.ilam(2).
%  v=upde normal unknowns part,   w=normal active uaux part
%  psi=u(p.nu/2+1:p.nu) pde part of adjoint kernel eigenvector
%  om=u(p.nu+p.naux+1:length(u)) active non-primary uaux part of kernel eigenvector 
%  G=normal pde equations, Q=normal auxiliary equations
%  Gv=partial derivative G_v, others analogous
%  Gvvpsi=\pa_v(G_v'*psi)
jtest=0; % switch for checking jac 
p=bpreduce(p); % set regular case sizes to compute G and Gu'*psii
nu=p.nu; nq=p.nc.nq; del=p.nc.del;
u1=[u(1:nu);u(2*nu+1:2*nu+p.naux)]; % u1(p.nu+1:end)' % pde-vars and pars 
nu1=length(u1); % org vars 
psipde=u(nu+1:2*nu); psiq=u(2*nu+p.naux+2:end); 
mu=u(2*nu+p.naux+1);  % 'zero' eval 
%r0=[r(1:nu); r(2*nu+1:2*nu+nq)]; % pde-resi and normal q
r0=resi(p,u1); % much better! 
[Gv,Glam]=getder(p,u1,r0); % get pde part linearization Glam=wrt old lam
if 0; jsw=p.sw.jac; p.sw.jac=0; Gv1=getGu(p,u1,r0); Jd=abs(Gv-Gv1); % test
    full(Jd(1:5,1:5)), p.sw.jac=jsw; pause; end
Gvpde=Gv(1:nu,1:nu); Glampde=Glam(1:p.nu); M=p.mat.M(1:nu,1:nu);
Gvpsi=Gvpde'*psipde; % reference value for Gu'*psi 
switch p.sw.spjac
    case 1; Gvvpsi=p.fuha.spjac(p,u);  % user def. 
    case 0; Gvvpsi=getGuTpsiu(p,u1,psipde);  % via numjac; BUT \pa_u(qu^T\psiti) still missing 
    otherwise; % pa_v(Gv*psi), expensive way 
  Gvvpsi=sparse(nu,nu);      
  for j=1:nu       
   up=u1-del*ej(j,nu1); r1=pderesi(p,up); Gv1=getGupde(p,up,r1); 
   r2=-(Gv1'*psipde-Gvpsi); r2=r2/del;      
   Gvvpsi=Gvvpsi+sparse(1:nu,j,r2,nu,nu); % add column with finite diff. approx.:   
   end   
end 
if jtest; Gvvpsi1=getGuTpsiu(p,u1,psipde); Jd=abs(Gvvpsi1-Gvvpsi); 
    full(Gvvpsi(1:4,1:4)), full(Gvvpsi1(1:4,1:4)),max(max(Jd)), pause; 
end 
% finite difference for mixed pde-aux part Gvwpsi=\pa_w(G_v*psi) and Gwvom=\pa_w(G_a*om)
uper=u1-del*ej(nu+p.nc.ilam(1),nu1); % perturb by old primary parameter
r2=pderesi(p,uper); % omitting mu*M*psipde; 
Gv2=getGupde(p,uper,r2);
Gvpsilam=-(Gv2-Gvpde)'*psipde/del;  % mixed deriv. 
if(nq==0)  
  % assemble Gu from block matrices; 3rd row is the normalization
  Gu=[[ Gv      mu*M     Glampde   M*psipde];  % G+mu*psi
      [Gvvpsi      Gv'   Gvpsilam  0*psipde]; % G^T*psi
      [0*psipde' 2*psipde'  0         0    ]     % |psi|^2-1
      [Gvpsilam'  Glam'     0         0    ]];   % Glam'*psi   
  return; 
end
% now case nq>0; exclude new primary parameter, and also old one,
% cause 'Gw' means derivative wrt OLD secondary one (HU) 
Gvpsiw=sparse(nu,nq); 
for k=1:nq 
   up=u1-del*ej(nu+p.nc.ilam(k+1),nu1); r1=pderesi(p,up); 
   Gv1=getGupde(p,up,r1); 
   r2=-(Gv1*psipde-Gvpsi)/del;
   Gvpsiw=Gvpsiw+sparse(1:nu,k,r2,nu,nq);
end
Gw=sparse(nu,nq); Qw=sparse(nq,nq); rq=r(2*nu+1:2*nu+nq); 
for k=1:nq    % Gw and Qw finite diff via pde:   
   k1=nu+p.nc.ilam(k+1); up=u1+del*ej(k1,nu1); 
   r1=pderesi(p,up); %+mu*M*psipde;    
   Gw=Gw+sparse(1:nu,k,(r1-r0(1:nu))/del,nu,nq);              
   r1=p.fuha.qf(p,up);
   Qw=Qw+sparse(1:nq,k,(r1-rq)/del,nq,nq);   
end
Glampde=Glam(1:nu); qlam=Glam(nu+1:end); 
zeroqmu=sparse(nq,1); 
zeroGw=sparse(nu,nq); zeroQw=sparse(nq,nq); 
qu=p.fuha.qfder(p,p.u); zero1q=sparse(1,nq); 
quupsi=p.fuha.quupsi(p,u1,r); % numerical quupsi still 'to do'
%  \pa_v       \pa_psi  pa_lam  \pa_w     \pa_mu    \pa_\psiq  % in H: 
z1=[Gvpde       mu*M   Glampde    Gw      M*psipde   zeroGw];  % G+mu*psi
z2=[Gvvpsi      Gvpde'  Gvpsilam Gvpsiw   0*psipde      qu'];  % Gu'*psi+qu'*psiq
z3=[qu           0*qu    qlam    Qw      zeroqmu     zeroQw];  % q+mu*psiq
z4=[quupsi       Gw'     0*qlam  0*Qw    zeroqmu        Qw'];  % Gw'*psi+qw'*psiq
z5=[0*psipde'  2*psipde' zero1q   0         0       2*psiq'];  % |psi|^2-1
z6=[Gvpsilam'  Glampde'  zero1q   0         0         qlam'];  % Glam'*psi+qlam'*psiq 
%size(z1),size(z2),size(z3),size(z4),size(z5),size(z6)
Gu=[z1;z2;z3;z4;z5;z6]; %mclf(6); spy(Gu); pause 
end 
