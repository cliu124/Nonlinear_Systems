% mod of diskpdeo to ellipse; r=radius, nr=#of boundary points, del=ellipticity
classdef diskpdeo< pde  
methods(Access = public)
  function o=diskpdeo(r,nr,del) % constructor 
     o.grid=grid2D; s=linspace(0,2*pi,nr); 
     o.grid.freeGeometry([r*cos(s)*del;r*sin(s)]); o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

