function [xg,c,ng]=discr_lob(xe,np)
% discr_lob:  mod of discr_lob from FSELIB _lo  
% xe: element end-nodes,  np=polynomial-order(1:np) 
% xg=all nodes, c=connectivity matrix, ng=number of global nodes
ne=size(xe,2)-1; 
for l=1:ne
  m=np(l);
  xi(1)=-1.0;
  if(m>1)
   [tL, wL]=lobatto(m-1);
   for j=2:m
    xi(j)=tL(j-1);
   end
  end
  xi(m+1)=1.0; mx=0.5*(xe(l+1)+xe(l)); dx=0.5*(xe(l+1)-xe(l));
  for j=1:m+1
    xen(l,j)=mx + xi(j)*dx;
  end
  for j=1:m+1
    xien(l,j)= xi(j);
  end
end
% connectivity matrix global nodes: 
Ic=2;      % global node counter
for l=1:ne
  Ic=Ic-1;
  for j=1:np(l)+1
    c(l,j)=Ic; xg(Ic)=xen(l,j);  Ic=Ic+1;
  end
end
ng=Ic-1;