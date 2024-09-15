function bjkl=getbjkl(p,u,phiv,del,nv)
% getbjkl: compute the directional derivatives bjkl=Guu[phi_j,n_kl],  m x m x m x n  matrix 
m=size(phiv,2); n=size(phiv,1); bjkl=zeros(m,m,m,n); deli=1/del; r=0*u(1:n); 
for j=1:m % compute Jacs at u+del*phiv(:,j)
  up=au2u(p,u2au(p,u)+del*phiv(:,j)); rp=resi(p,up); 
  um=au2u(p,u2au(p,u)-del*phiv(:,j)); rm=resi(p,um); 
  Gup=getGu(p,up,rp); Gum=getGu(p,um,rm); Gud=0.5*deli*(Gup-Gum); 
  for k=1:m
     for l=1:m
         bjkl(j,k,l,:)=Gud*reshape(nv(k,l,:),n,1); 
     end
  end
end 