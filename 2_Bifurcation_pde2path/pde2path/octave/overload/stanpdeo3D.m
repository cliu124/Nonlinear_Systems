classdef stanpdeo3D < pde 
% stanpdeo3D:  standard 3D OOPDE object for p2p (classdef), adapted to OCTAVE 
% nargin=4: lx,ly,lz, h;  nargin=5: lx,ly,lz,h, sw; 
% nargin=6: lx,ly,lz, nx,ny,nz;  nargin=7: lx,ly,lz, nx,ny,nz, sw; 
% sw.ref: #ref.steps, sw.sym: if 1, then use ccbar ("criss-cross" like mesh) 
methods(Access = public)
  function o=stanpdeo3D(lx,ly,lz,varargin) % constructor 
     o.grid=grid3D;  h=varargin{1}; sw=[]; 
     if nargin==5; sw=varargin{2}; end; if nargin==7; sw=varargin{4}; end
     try sym=sw.sym; catch; sym=0; end;  
     sym 
     switch sym 
      case -1; o.grid.myCube(lx,ly,lz,h); % the simplest case 
      case 0; o.grid.bar(-lx,lx,-ly,ly,-lz,lz,h); 
      case 1; o.grid.ccbar(-lx,lx,-ly,ly,-lz,lz,h); 
     end 
     try ref=sw.ref; catch; ref=0; end 
     for i=1:ref; o.grid.refineMesh; end;  
     o.fem=lagrange13D; 
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
