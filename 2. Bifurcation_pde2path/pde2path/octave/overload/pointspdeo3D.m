classdef pointspdeo3D < pde 
% stanpdeo3D:  standard 3D OOPDE object for p2p (classdef), adapted to OCTAVE 
% nargin=4: lx,ly,lz, h;  nargin=5: lx,ly,lz,h, sw; 
% nargin=6: lx,ly,lz, nx,ny,nz;  nargin=7: lx,ly,lz, nx,ny,nz, sw; 
% sw.ref: #ref.steps, sw.sym: if 1, then use ccbar ("criss-cross" like mesh) 
methods(Access = public)
  function o=pointspdeo3D(po,varargin) % constructor 
     o.grid=grid3D;   
     o.grid.p2dom(po); 
     try ref=sw.ref; catch; ref=0; end 
     for i=1:ref; o.grid.refineMesh; end;  
     o.fem=lagrange13D; 
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
