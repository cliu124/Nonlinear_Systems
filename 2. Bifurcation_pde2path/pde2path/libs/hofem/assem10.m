function [K,M]=assem10(p)
% assem10: assemble K and M for 10-nodes tetras 
% Adapted from Pozrikidis' FSElib 
% not vectorized, slow. ng=#nodes, p=point coordinates, 
% c=nelem x 10 connectivity matrix 
c=p.hofem.tri; po=p.pdeo.grid.p'; ng=size(po,1); % read 10-node tetras, and points from p
gdm=0*speye(ng); gmm=gdm; % initialize
ne=size(c,1); volume=0; %ng, ne, size(gdm), pause 
for l=1:ne     % loop over the elements
    % compute the element diffusion matrix
j=c(l,1); x1=po(j,1); y1=po(j,2); z1=po(j,3);
j=c(l,2); x2=po(j,1); y2=po(j,2); z2=po(j,3);
j=c(l,3); x3=po(j,1); y3=po(j,2); z3=po(j,3);
j=c(l,4); x4=po(j,1); y4=po(j,2); z4=po(j,3);
j=c(l,5); x5=po(j,1); y5=po(j,2); z5=po(j,3);
j=c(l,6); x6=po(j,1); y6=po(j,2); z6=po(j,3);
j=c(l,7); x7=po(j,1); y7=po(j,2); z7=po(j,3);
j=c(l,8); x8=po(j,1); y8=po(j,2); z8=po(j,3);
j=c(l,9); x9=po(j,1); y9=po(j,2); z9=po(j,3);
j=c(l,10);x10= po(j,1); y10= po(j,2); z10= po(j,3);

[edm_elm, emm_elm, vlm_elm]=edmm_t10(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, x5,y5,z5, x6,y6,z6 ...
   ,x7,y7,z7, x8,y8,z8, x9,y9,z9, x10,y10,z10);
 volume=volume + vlm_elm;
  for i=1:10
     i1=c(l,i);
     for j=1:10
       j1=c(l,j);
       gdm(i1,j1)=gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1) = gmm(i1,j1) + emm_elm(i,j);
     end
  end
end
K=sparse(gdm); M=sparse(gmm); 
