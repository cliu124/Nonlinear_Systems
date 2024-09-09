function [fij,gij]=getfg(p,u,phiv,psiv,chi,del) 
% getfg: quadratic coefficients in bif.equations 
% f_ij=psi_i*G_{u,lam}*phi_j, g_ij=psi_i*Guu[phi_j,chi] 
m=size(phiv,2); fij=zeros(m,m); gij=zeros(m,m); deli=1/del; 
for i=1:m
  for j=1:m    
     up=au2u(p,u2au(p,u)+del*phiv(:,j)); rp=resi(p,up); 
     um=au2u(p,u2au(p,u)-del*phiv(:,j)); rm=resi(p,um); 
     [Gup,Glamp]=getder(p,up,rp); [Gum,Glamm]=getder(p,um,rm); 
     Glamd=0.5*deli*(Glamp-Glamm);  Gud=0.5*deli*(Gup-Gum); % central diff
     fij(i,j)=psiv(:,i)'*Glamd;  gij(i,j)=psiv(:,i)'*(Gud*chi); 
  end
end 