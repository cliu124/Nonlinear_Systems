classdef stanpdeo3D < pde 
% stanpdeo3D:  standard 3D OOPDE object for p2p (classdef) 
% nargin=4: lx,ly,lz, h;  nargin=5: lx,ly,lz,h, sw; 
% nargin=6: lx,ly,lz, nx,ny,nz;  nargin=7: lx,ly,lz, nx,ny,nz, sw; 
% sw.ref: #ref.steps, sw.sym: if 1, then use ccbar ("criss-cross" like mesh) 
methods(Access = public)
  function o=stanpdeo3D(lx,ly,lz,varargin) % constructor 
     o.grid=grid3D;  ref=0; sym=0;
     if nargin==5; sw=varargin{2}; end; if nargin==7; sw=varargin{4}; end
     try sym=sw.sym; catch; end;  try ref=sw.ref; catch; end 
     if nargin==4 || nargin==5;          
        if ~sym  o.grid.bar(-lx,lx,-ly,ly,-lz,lz,varargin{1});  
        else  o.grid.ccbar(-lx,lx,-ly,ly,-lz,lz,varargin{1});
        end
     else
        if ~sym  o.grid.bar(-lx,lx,-ly,ly,-lz,lz,varargin{1},varargin{2},varargin{3});
        else  o.grid.ccbar(-lx,lx,-ly,ly,-lz,lz,varargin{1},varargin{2},varargin{3});
        end
     end
     for i=1:ref; o.grid.refineMesh; end; % can make the mesh more symmetric
     o.fem=lagrange13D; 
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
