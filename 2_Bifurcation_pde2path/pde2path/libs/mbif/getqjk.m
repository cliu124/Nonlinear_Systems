function qjk=getqjk(p,Gu,u,phiv,del)
% getqjk: compute the directional derivatives qjk=Guu[phi_j,phi_k] 
% m x m matrix of n vectors 
m=size(phiv,2); n=size(phiv,1); qjk=zeros(m,m,n); deli=1/del; 
for j=1:m
    up=au2u(p,u2au(p,u)+del*phiv(:,j)); rp=resi(p,up); Gup=getGu(p,up,rp); 
    um=au2u(p,u2au(p,u)-del*phiv(:,j)); rm=resi(p,um); Gum=getGu(p,um,rm); 
    Gud=0.5*deli*(Gup-Gum); % central differences
    for k=1:m; qjk(j,k,:)=Gud*phiv(:,k); end
end

