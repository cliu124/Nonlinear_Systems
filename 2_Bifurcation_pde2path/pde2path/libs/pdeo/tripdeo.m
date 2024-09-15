classdef tripdeo< pde  
% tripdeo (classdef): domain and mesh for an equilat. triangle 
%  pdeo=tripdeo(r,hmax)

methods(Access = public)
  function o=tripdeo(r,hmax) % constructor 
     o.grid=grid2D; s3=sqrt(3); 
    % x=[-0.5 0.5 0]; y=[0 0 s3/2];  
     x=[-0.5 0 0.5 0.25 0 -0.25]; y=[0 0 0 s3/4 s3/2 s3/4];  
     o.grid.freeGeometry(r*[x;y])
     h=1; while h > hmax; o.grid.refineMesh; h = h/2; end
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

