function L=cotmatrix(V,F)
  % COTMATRIX computes cotangent matrix (laplacian mesh operator), (mass/area
  % terms already cancelled out)
  %
  % L=cotmatrix(V,F)
  % L=cotmatrix(V,T)
  %
  % For size(F,2)==4, This is distinctly NOT following definition that appears
  % in the appendix of: ``Interactive Topology-aware Surface Reconstruction,''
  % by Sharf, A. et al
  % http://www.cs.bgu.ac.il/~asharf/Projects/InSuRe/Insure_siggraph_final.pdf
  %
  % Instead it is a purely geometric construction. Find more details in Section
  % 1.1 of "Algorithms and Interfaces for Real-Time Deformation of 2D and 3D
  % shapes" [Jacobson 2013]
  %
  % ND derivation given in "A MONOTONE FINITE ELEMENT SCHEME FOR
  % CONVECTION-DIFFUSION EQUATIONS" [Xu & ZIKATANOV 1999]
  %
  % 3D derivation given in "Aspects of unstructured grids and finite-volume
  % solvers for the Euler and Navier-Stokes equations" [Barth 1992]
  %
  %
  % Inputs:
  %   V  #V x dim matrix of vertex coordinates
  %   F  #F x simplex-size matrix of indices of triangle or tetrahedron corners
  % Outputs:
  %   L  sparse #V x #V matrix of cot weights 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Denis Zorin
  %
  % See also: cotangent
  %

  ss=size(F,2);
  switch ss
  case 2
    G=grad(V,F);
    A=diag(sparse(edge_lengths(V,F)));
    L=-G'*A*G;
  case 3
    %% Could just replace everything with:
    if 1
    C=cotangent(V,F);
    L=sparse(F(:,[2 3 1]), F(:,[3 1 2]), C,size(V,1),size(V,1));
    L=L+L';
    L=L-diag(sum(L,2));
    else 
    if(size(F,1) == 3); warning('F seems to be 3 by #F, it should be #F by 3'); end
    F=F'; i1=F(1,:); i2=F(2,:); i3=F(3,:); 
    % #F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1=V(i3,:) - V(i2,:);  v2=V(i1,:) - V(i3,:); v3=V(i2,:) - V(i1,:);
    % computing *unsigned* areas 
    if size(V,2) == 2
        dblA=abs(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1));
    elseif size(V,2) == 3
        %n =cross(v1,v2,2);  dblA =multinorm(n,2);
        % area of parallelogram is twice area of triangle       
        n =cross(v1,v2,2); 
        dblA=(sqrt(sum((n').^2)))';
    else 
        error('unsupported vertex dimension %d', size(V,2))
    end
    % cotangents and diagonal entries for element matrices
    cot12=-dot(v1,v2,2)./dblA/2; cot23=-dot(v2,v3,2)./dblA/2; cot31=-dot(v3,v1,2)./dblA/2;
    % diag entries computed from the condition that rows of the matrix sum up to 1
    % (follows from  the element matrix formula E_{ij}=(v_i dot v_j)/4/A )
    diag1=-cot12-cot31; diag2=-cot12-cot23; diag3=-cot31-cot23;
    % indices of nonzero elements in the matrix for sparse() constructor
    i=[i1 i2 i2 i3 i3 i1  i1 i2 i3];
    j=[i2 i1 i3 i2 i1 i3  i1 i2 i3];
    % values corresponding to pairs form (i,j)
    v=[cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
    % for repeated indices (i,j) sparse automatically sums up elements, as we
    % want
    L=sparse(i,j,v,size(V,1),size(V,1));
    end
  end
end
