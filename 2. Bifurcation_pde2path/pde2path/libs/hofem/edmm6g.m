function [edm,emm,arel]=edmm6g(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,NQ,c,isw)
% edmm6g: diffusion and mass matrices for a 6-node triangle (FSELib) 
% general c=[c1(:) c2(:); c3(:) c4(:)],1:6,
% c1=\pa_x^2,c4=\pa_y^2,c2,3 mixed
% compute the mapping coefficients
[al,be,ga]=elm6_abc(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6);
[xi,eta,w]=gauss_trgl(NQ); % read the triangle quadrature
[n,m]=size(c); 
if n==1; % c=one row
   if m==1; c1=c*ones(1,6); c4=c1; % scalar 
   else c1=c(1:6); c4=c1; end 
else % 2 rows 
   if m==1; c1=c(1,1)*ones(1,6); c4=c(2,2)*ones(1,6);
   else c1=c(1,1:6); c2=c(1,7:12); c3=c(2,1:6); c4=c(2,7:12); end 
end 
edm=zeros(6); emm=edm; 
% perform the quadrature
arel=0.0;  % element area (optional)
for i=1:NQ
 [psi,gpsi,hs]=elm6_interp(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,al,be,ga,xi(i),eta(i));
 cf=0.5*hs*w(i); %size(c1), size(psi), pause 
 c1i=dot(c1,psi); c2i=dot(c2,psi); c3i=dot(c3,psi); c4i=dot(c4,psi); 
 for k=1:6
  for l=1:6
if ~isw % use row wise c 
   edm(k,l)=edm(k,l)+cf*c1(k)*gpsi(k,1)*gpsi(l,1)+cf*c4(k)*gpsi(k,2)*gpsi(l,2) ...
        +cf*c2(k)*gpsi(k,2)*gpsi(l,1)+cf*c3(k)*gpsi(k,1)*gpsi(l,2);  
else % use interpol. c
    edm(k,l)=edm(k,l)+cf*c1i*gpsi(k,1)*gpsi(l,1)+cf*c4i*gpsi(k,2)*gpsi(l,2) ...
             +cf*c2i*gpsi(k,2)*gpsi(l,1)+cf*c3i*gpsi(k,1)*gpsi(l,2); 
end
    emm(k,l)=emm(k,l)+psi(k)*psi(l)*cf;
  end
 end
 arel=arel+cf;
end