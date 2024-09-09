classdef hexpdeo< pde  
% hexpdeo (classdef): domain and mesh for a unit hexagon 
methods(Access = public)
  function o=hexpdeo(hmax) % constructor 
     o.grid=grid2D; s3=sqrt(3); 
     x=[1 0.5 -0.5 -1 -0.5 0.5]; y=[0 s3/2 s3/2 0 -s3/2 -s3/2];  
     o.grid.freeGeometry([x;y])
     h=1; while h > hmax; o.grid.refineMesh; h = h/2; end
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

