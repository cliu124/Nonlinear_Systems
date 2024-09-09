classdef wankelpdeo< pde  
% wankelpdeo: domain and mesh for an equilat.triangle with rounded corners
methods(Access = public)
  function o=wankelpdeo(nbd,hmax) % constructor 
     o.grid=grid2D; t=linspace(0,2*pi,nbd); t=t(1:end-1); R=1; r=R/6*2; a=0.5;  
    x=(R-r)*cos(t)+a*r*cos((R/r-1)*t); 
    y=(R-r)*sin(t)-a*r*sin((R/r-1)*t); 
     o.grid.freeGeometry([x;y])
     h=1; while h > hmax; o.grid.refineMesh; h = h/2; end
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

