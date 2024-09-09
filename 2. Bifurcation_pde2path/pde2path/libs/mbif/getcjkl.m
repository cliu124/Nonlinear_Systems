function cjkl=getcjkl(p,u,phiv,del)
% getcjkl: compute the directional derivatives cjkl=Guuu[phi_j,phi_j,phi_k]
% m x m x m matrix of n vectors 
m=size(phiv,2); n=size(phiv,1); cjkl=zeros(m,m,m,n); deli2=1/del^2; 
for j=1:m
  for k=1:m
    upp=au2u(p,u2au(p,u)+del*(phiv(:,j)+phiv(:,k))); 
    upm=au2u(p,u2au(p,u)+del*(phiv(:,j)-phiv(:,k))); 
    ump=au2u(p,u2au(p,u)+del*(-phiv(:,j)+phiv(:,k))); 
    umm=au2u(p,u2au(p,u)-del*(phiv(:,j)+phiv(:,k))); 
    rpp=resi(p,upp);  rpm=resi(p,upm); 
    rmp=resi(p,ump);  rmm=resi(p,umm); 
    Gupp=getGu(p,upp,rpp); Gupm=getGu(p,upm,rpm);  
    Gump=getGu(p,ump,rmp); Gumm=getGu(p,umm,rmm);  
    Gud=Gupp-Gupm-Gump+Gumm; 
    %deli2*max(max(abs(Gud))), pause 
    for l=1:m; cjkl(j,k,l,:)=0.25*deli2*(Gud*phiv(:,l)); end
  end
end

