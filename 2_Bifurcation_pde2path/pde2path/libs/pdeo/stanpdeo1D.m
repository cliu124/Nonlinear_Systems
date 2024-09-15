classdef stanpdeo1D < pde 
% stanpdeo1D:  standard 1D OOPDE object for p2p (classdef) 
%
% pdeo=stanpdeo1D(lx,h)

methods(Access = public)
  function o=stanpdeo1D(lx,h) % constructor 
      o.grid=grid1D; o.grid.interval([-lx,lx],h); o.fem=lagrange11D;
  end
end  
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
