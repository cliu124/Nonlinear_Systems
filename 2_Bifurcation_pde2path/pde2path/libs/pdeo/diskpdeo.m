classdef diskpdeo< pde  
% diskpdeo: disk pde object, arguments R and nphi (#of points on circle) 
%
%  pdeo=diskpdeo(R,npphi)
% R=radius, nphi=#points on radius 

methods(Access = public)
  function o=diskpdeo(R,nphi) % constructor 
     o.grid=grid2D; s=linspace(0,2*pi,nphi);  
     o.grid.freeGeometry([R*cos(s);R*sin(s)]); 
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

