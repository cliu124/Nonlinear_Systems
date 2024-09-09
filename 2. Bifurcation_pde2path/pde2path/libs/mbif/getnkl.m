function [nv,chi,p]=getnkl(p,u,Gu,Glam,phiv,psiv,del) 
% getnkl: project quadratic terms in bif.eqns back onto kernel 
%  nv(k,l)=-0.5*P(Gu\Guu[phi_k,phi_l]); chi=-P(Gu\Glam)
m=size(phiv,2); n=size(phiv,1); qjk=getqjk(p,Gu,u,phiv,del); nv=zeros(m,m,n); 
tic
for k=1:m 
  for l=1:m
    [nt,p]=p.fuha.lss(Gu,reshape(qjk(k,l,:),n,1),p); 
     nv(k,l,:)=-0.5*proNperp(nt,psiv,phiv); 
  end
end
toc
chi=Gu\Glam; chi=-proNperp(chi,psiv,phiv); 