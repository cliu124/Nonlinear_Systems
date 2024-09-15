function pot=pot(p) % potential 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; pot=x.^2+y.^2;
