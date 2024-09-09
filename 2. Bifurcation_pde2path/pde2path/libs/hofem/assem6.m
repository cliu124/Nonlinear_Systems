function [K,M]=assem6(p)
% assem6: assemble K and M for 6nodes triangles, here constant coeff, 
% adapted from Pozrikidis' FSElib 
c=p.hofem.tri; po=p.pdeo.grid.p';  % read 6-node triang. and points from p
ng=size(po,1); gdm=0*speye(ng); gmm=gdm; % initialize
NQ=12; ne=size(c,1); 
for l=1:ne     % loop over the elements
% compute the element diffusion and mass matrices
j=c(l,1); x1=po(j,1); y1=po(j,2);
j=c(l,2); x2=po(j,1); y2=po(j,2);
j=c(l,3); x3=po(j,1); y3=po(j,2);
j=c(l,4); x4=po(j,1); y4=po(j,2);
j=c(l,5); x5=po(j,1); y5=po(j,2);
j=c(l,6); x6=po(j,1); y6=po(j,2);

[edm_elm, emm_elm, arel]=edmm6(x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ);

  for i=1:6
     i1=c(l,i);
     for j=1:6
       j1=c(l,j);
       gdm(i1,j1)=gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1)=gmm(i1,j1) + emm_elm(i,j);
     end
   end
end
K=sparse(gdm); M=sparse(gmm); 
