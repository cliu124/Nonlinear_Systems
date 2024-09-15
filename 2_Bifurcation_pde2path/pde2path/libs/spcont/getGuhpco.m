function Hu=getGuhpco(p,u,r)
% getGufoco: HP continuation Jacobian 
%
% see hpcresi for the pertinent extended system H 
%     hpcontini for the setup of ilam
%     hpreduce for reducing to the original PDE+constraints
%
% build H_u from blocks 
% p.sw.spjac=1  highly recommended for speed 
% some second mixed derivatives ignored here (quwr, quwi,.. below) 
% because usually all zero. See getGufoco and put in if needed!
p=hpreduce(p); % set regular case sizes to compute G and phi
nu=p.nu; nq=p.nc.nq; 
M=p.mat.M(1:nu,1:nu); del=p.nc.del; ilam=p.nc.ilam; 
u1=[u(1:nu); u(3*nu+1:3*nu+p.naux)]; % pde & pars  (no phi,om)
nu1=length(u1);  % org variables (no phi,om)
phir=u(nu+1:2*nu); phii=u(2*nu+1:3*nu); om=u(3*nu+p.naux+1); % PDE-evecs, and om 
r1=r(1:nu); rq=r(nu+1:nu+nq); % resis
%r1=p.fuha.sG(p,u1); %rq=p.fuha.qf(p,u1); % alternative
Gv=getGupde(p,u1,r1); 
Gvr=Gv*phir; Gvi=Gv*phii; 
switch p.sw.spjac
    case 1; Gvvph=p.fuha.spjac(p,u); % user def. 
    case 0; Gvvph=getGuphih(p,u1,phir,phii); % numjac 
    otherwise ; Gvvph=sparse(2*nu,nu); % expensive way
 for j=1:nu  
  up=u1-del*ej(j,nu1); rp=pderesi(p,up); % backward 
  Gvp=getGupde(p,up,rp); % perturbed pde-part lin.
  rr=-(Gvp*phir-Gvr)/del; ri=-(Gvp*phii-Gvi)/del;   
  Gvvph=Gvvph+sparse(1:2*nu,j,[rr;ri],2*nu,nu); 
 end  
end
% perturb by OLD primary parameter lam (other colums (nq>0) later) 
up=u1; up(nu+ilam(1))=up(nu+ilam(1))+del; 
rp=pderesi(p,up); Gv1=getGupde(p,up,rp); 
Glam=(rp-r1)/del;  
H2lam=(Gv1*phir-Gvr)/del; 
H3lam=(Gv1*phii-Gvi)/del; 
if(nq==0)  % the simple case
 %     pa_u            pa_phir pa_phii pa_w    pa_om
  Hu=[[ Gv                0*Gv   0*Gv   Glam     0*phir]; 
      [Gvvph(1:nu,:)      Gv     om*M   H2lam   M*phii]; 
      [Gvvph(nu+1:2*nu,:) -om*M  Gv     H3lam    -M*phir]; 
      [0*phir'             p.c  0*p.c    0      0  ];
      [0*phir'            0*p.c  p.c     0      0  ]]; 
  return
end
% case nq>0: keep Glam, H2lam, H3lam, further OLD aux pars at p.nc.ilam(2:nq+1), 
% der.wrt these called Gw, H2w, etc 
% quwr=\pa_w(qu*phir), quwi=\pa_w(qu*phii)  etc 
%Gv=Gv(1:nu,1:nu);  % PDE part
Gw=sparse(nu,nq); H2w=Gw; H3w=Gw; qw=sparse(nq,nq); 
qulamr=sparse(nq,1); qulami=qulamr; quwr=sparse(nq,nq); quwi=quwr; % still to do! 
%p.nc.ilam(1:nq+1), pause 
%H5=r(3*nu+nq+1:3*nu+2*nq); H6=r(3*nu+2*nq+1:3*nu+3*nq);
rp0=p.fuha.sG(p,u1); % recompute pderesi, (better than taking from r) 
rq0=p.fuha.qf(p,u1); % norm(rq0-rq), pause 
upert=u1+p.nc.del*ej(nu+p.nc.ilam(1),nu1); 
qlam=(p.fuha.qf(p,upert)-rq)/del; 
for i=1:nq % derivatives of G and q wrt old aux. vars 
    upert=u1+p.nc.del*ej(nu+p.nc.ilam(i+1),nu1); 
    rp1=p.fuha.sG(p,upert); 
    rq1=p.fuha.qf(p,upert);
    Gw(:,i)=(rp1-rp0)/del; 
    qw(:,i)=(rq1-rq0)/del;
    Gup=getGupde(p,upert,[rp1;rq1]); 
    if 0; jsw=p.sw.jac; p.sw.jac=1; Gup2=getGupde(p,upert,[rp1;rq1]); 
        Gud1=abs(Gup-Gup2); max(max(Gud1)), pause; p.sw.jac=jsw; end 
    Gud=Gup-Gv; 
    H2w(:,i)=Gud*phir; 
    H3w(:,i)=Gud*phii;         
end
cv=p.c(1:nu); cvt=p.c(nu+1:end); 
phir=phir(1:nu); phii=phii(1:nu);  % pde-part of Evec 
zuu=sparse(nu,nu); zuq=sparse(nu,nq); % zero nu x nu and nu x nq
quupr=p.fuha.quuphir(p,u,r); quupi=p.fuha.quuphii(p,u,r); 
zqq=sparse(nq,nq); % zero nq x nq
zq1=sparse(nq,1); 
qu=p.fuha.qfder(p,up); 
%   pa_u           pa_phir  pa_phii  pa_lam pa_w   pa_om  pa_phitr pa_phiti
z1=[Gv                zuu    zuu     Glam     Gw    0*phir     zuq  zuq]; 
z2=[Gvvph(1:nu,:)      Gv     om*M   H2lam    H2w    M*phii    Gw    zuq]; 
z3=[Gvvph(nu+1:2*nu,:) -om*M  Gv     H3lam    H3w   -M*phir    zuq   Gw]; 
z4=[qu                 zuq'   zuq'   qlam     qw      zq1       zqq   zqq]; 
z5=[quupr             qu     zuq'    qulamr   quwr   zq1       qw   zqq]; 
z6=[quupi             zuq'   qu     qulami   quwi   zq1       zqq   qw]; 
z7=[zuq'               cv    0*cv      0      zq1'   zq1'      cvt  zq1'];
z8=[zuq'               0*cv   cv       0      zq1'   zq1'      zq1'  cvt]; 
Hu=[z1; z2; z3; z4; z5; z6; z7; z8]; %mclf(6); spy(Gu); size(Gu), pause 