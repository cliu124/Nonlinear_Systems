function r=hpcresi(p,u)
% hpcresi: residual for extended system for HP continuation
%
% H(U)=[G; Gu*phir+om*M*phii; Gu*phi-om*M*phir, c'*phir-1; c'*phii]
% 
% dim(H)=3(nu+nq)+
% 
% ||phi||^2-1]; \in R^(2*nu+2*nq+1) 
% 
% or rather 
%
% 
% where: 
% G=regular problem part=[Gpde;q] (dimension nu+nq)
% phi=phir+i*phii=kernel eigenvector  (dimension nu+nq)
% Gu=derivative of G (incl.q) wrt (u,uaux)  (dimension nu+nq x nu+nq) 
% 
p=hpreduce(p); % turn off hpcont and remove new foco primary param 
nu=p.nu; nq=p.nc.nq; 
u1=[u(1:nu); u(3*nu+1:3*nu+1+p.naux)]; nu1=length(u1); % normal pde part and pars 
rp=pderesi(p,u1); if(nq>0) rq=p.fuha.qf(p,u1); else rq=[]; end % pde and aux part:
r1=[rp;rq]; cv=p.c; %(1:nu); 
% append further resis: 
Gu=getGu(p,u1,r1); Gupde=Gu(1:nu,1:nu); 
M=p.mat.M(1:nu,1:nu); 
phir=u(nu+1:2*nu); phii=u(2*nu+1:3*nu); 
om=u(3*nu+p.naux+1); 
r2=Gupde*phir+om*M*phii; r3=Gupde*phii-om*M*phir; 
if nq==0; % sort into residual, simple case 
    r4=cv*phir-1; r5=cv*phii; 
    r=[rp; r2; r3; r4; r5]; return; 
end 
% case nq>0, more complicated 
qu=p.fuha.qfder(p,u1); 
Gw=sparse(nu,nq); qw=sparse(nq,nq); %p.nc.ilam,nq
for i=1:nq % derivatives of G and q wrt old aux. vars (incl.old prim)
    upert=u1+p.nc.del*ej(p.nu+p.nc.ilam(i+1),nu1); 
    rp1=p.fuha.sG(p,upert); 
    rq1=p.fuha.qf(p,upert);
    Gw(:,i)=(rp1-rp)/p.nc.del; 
    qw(:,i)=(rq1-rq)/p.nc.del;
end
phits=3*nu+p.naux+2; % startindex of q-part of Evec 
phirt=u(phits:phits+nq-1); phiit=u(phits+nq:end); % q-part(s) of Evec 
r2=r2+Gw*phirt; r3=r3+Gw*phiit; 
r5=qu*phir+qw*phirt; r6=qu*phii+qw*phiit;
r7=cv*[phir;phirt]-1; r8=cv*[phii;phiit];  
r=[rp; r2; r3; rq; r5; r6; r7; r8]; 
return % uncomment to get infos
qw, cv(end-5:end), phirt, phiit
rv=[norm(rp,'inf'), norm(r2,'inf'), norm(r3,'inf'), norm(rq,'inf'), norm(r5,'inf'), ...
    norm(r6,'inf'), norm(r7,'inf'),  norm(r8,'inf')]; rv, pause 
