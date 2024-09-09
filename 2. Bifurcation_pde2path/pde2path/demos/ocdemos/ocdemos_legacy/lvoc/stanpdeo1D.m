%% 
% stanpdeo1D:  standard 1D OOPDE object for p2p (classdef) 
classdef stanpdeo1D < pde 
methods(Access = public)
  function o=stanpdeo1D(lx,h) % constructor 
      o.grid=grid1D; o.grid.interval([0,lx],h); o.fem=lagrange11D;
  end
end  
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
