classdef freegeompdeo< pde  
% freegeompdeo: free geometry, x,y=bd-coordinates, hmax=max mesh-width
%   not caring for boundary segment identifiers (all 1) 
methods(Access = public)
  function o=freegeompdeo(x,y,hmax) % constructor 
     o.grid=grid2D;  o.grid.freeGeometry([x;y]); 
     h=1; while h > hmax; o.grid.refineMesh; h = h/2; end
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

