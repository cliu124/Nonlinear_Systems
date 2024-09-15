classdef stanpdeo1Db < pde 
% stanpdeo1Db: 1D pdeo, x\in(l1,l2) (not x\in(-lx,lx) as in stanpdeo1D) 
methods(Access = public)
  function o=stanpdeo1Db(l1,l2,h) % constructor 
      o.grid=grid1D; o.grid.interval([l1,l2],h); o.fem=lagrange11D;
  end
end  
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
