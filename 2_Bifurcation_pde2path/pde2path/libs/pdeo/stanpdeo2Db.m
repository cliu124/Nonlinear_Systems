classdef stanpdeo2Db < pde  
% stanpdeo2Db:  standard 2D OOPDE object for p2p (classdef) 
% nargin=3: lx,ly, h;   
% nargin=4: lx,ly, nx, ny; or lx, ly, h, sw;  
% nargin=5: lx,ly, nx, ny, sw. 
% sw.ref: #ref.steps, sw.sym: if 1, then use ccsquare (like sympoi)  
methods(Access = public)
  function o=stanpdeo2Db(lx,ly,varargin) % constructor 
     o.grid=grid2D; ref=0; sym=0;
     if nargin==3; o.grid.square(0,lx,0,ly,varargin{1}); end 
     if nargin==4; % lx,ly, nx, ny; or lx, ly, h, sw;
        if ~isstruct(varargin{2}); o.grid.square(0,lx,0,ly,varargin{1},varargin{2});
        else hmax=varargin{1}; sw=varargin{2}; 
            try sym=sw.sym; catch; end 
            switch sym; % 
                case 0; o.grid.square(0,lx,0,ly,hmax); 
                case 1; o.grid.ccsquare(0,lx,0,ly,hmax); 
                case 2;  % true criss-cross via sympoi
                    o.grid.square(0,lx,0,ly,hmax); 
                     rlongs=o.grid.rlong;  o.grid.rlong=1;
                idx=1:size(o.grid.t,2); % to refine all 
                 o.grid.refineMesh(idx); o.grid.rlong=rlongs;                              
            end
        end
     end 
     if nargin==5;  % lx,ly, nx, ny, sw. 
        nx=varargin{1}; ny=varargin{2}; sw=varargin{3}; 
        try sym=sw.sym; catch; end
         switch sym; % 
            case 0; o.grid.square(0,lx,0,ly,nx,ny); 
            case 1; o.grid.ccsquare(0,lx,0,ly,nx,ny); 
            case 2; p.plot.ifig=6; p.nc.neq=1; % fill p with some dummys  
                o.grid.square(0,lx,0,ly,nx,ny); 
                rlongs=o.grid.rlong;  o.grid.rlong=1;
                idx=1:size(o.grid.t,2); % to refine all 
                 o.grid.refineMesh(idx); o.grid.rlong=rlongs;                
         end
     end 
     try ref=sw.ref; catch; end  
     for i=1:ref; o.grid.refineMesh; end;
     o.fem=lagrange12D;
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
