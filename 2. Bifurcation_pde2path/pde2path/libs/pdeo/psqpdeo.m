classdef psqpdeo< pde  
% psqpde: domain and mesh for a perturbed square (rounded corners) 
methods(Access = public)
  function o=psqpdeo(nbd,hmax,del) % constructor 
     o.grid=grid2D; t=linspace(-pi,pi,nbd); t=t(1:end-1); 
     x1=t; y1=-pi-del*cos(t/2); x3=-t; y3=pi+del*cos(t/2); 
     x2=pi+del*cos(t/2); y2=t;  x4=-pi-del*cos(t/2); y4=-t; 
     x=[x1,x2,x3,x4]; y=[y1,y2,y3,y4];
     o.grid.freeGeometry([x;y]); 
     h=1; while h > hmax; o.grid.refineMesh; h = h/2; end
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

