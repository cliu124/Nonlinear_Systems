classdef stanpdeo2D < pde  
% stanpdeo2D:  standard 2D OOPDE object for p2p (classdef), adapted to OCTAVE
% nargin=3: lx,ly, h;   
% nargin=4: lx,ly, nx, ny; or lx, ly, h, sw;  
% nargin=5: lx,ly, nx, ny, sw. 
% sw.ref: #ref.steps, sw.sym: if 1, then use ccsquare (like sympoi)  
methods(Access = public)
  function o=stanpdeo2D(lx,ly,varargin) % constructor 
     o.grid=grid2D; try; sw=varargin{2};  catch; sw=[]; end; 
     try sym=sw.sym; catch; sym=0; end 
     try; ref=sw.ref;  catch; ref=0; end; 
    o.grid.square(-lx,lx,-ly,ly,varargin{1});      
   %  o.grid.mySquare(lx,ly,varargin{1}); 
   rlongs=o.grid.rlong;
     if sym==2;  o.grid.rlong=1;  end 
     if ref>0; 
       for i=1:ref; idx=1:size(o.grid.t,2);  o.grid.refineMesh(idx);   end 
     end 
     o.grid.rlong=rlongs; 
     o.fem=lagrange12D; 
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
