classdef secpdeo< pde  
% secpdeo: sector pde object, with radial mesh, based on freeGeoPts
%
% pdeo=secpdeo1(r,phi,nr,varargin) 
%
% r=radius, phi=half-angle, nr=#point in r 
% varargin.al=scaling of r-distr.of points
% varargin.phifac=factor for # of angular points
methods(Access = public)
  function o=secpdeo(R,phi,nr,varargin) % constructor 
     if nargin>3; al=varargin{1}; else al=1; end 
     if nargin>4; phifac=varargin{2}; else phifac=round(phi*3); end 
     o.grid=grid2DHU;  
     bda=linspace(0,R.^(1/al),nr); bdb=R.^(1/al)-bda(1:end-1); 
     nphi=(nr+1)*phifac+1; s=linspace(-phi,phi,nphi); s=s(2:end-1); 
     bdsi=[nr, nr+nphi-1]; % the splitting points for bd-segment identific. 
     bd=[bda.^al*cos(-phi) R*cos(s) cos(phi)*bdb.^al; 
         bda.^al*sin(-phi) R*sin(s) sin(phi)*bdb.^al]; 
     rv=bda.^al; x=[]; y=[]; % the inner mesh points, generated radially symmetric 
     for i=2:nr-1;  
        j=round(phifac*i)+1; 
        th=linspace(-phi,phi,j); th=th(2:end-1);     
        x=[x rv(i)*cos(th)]; y=[y rv(i)*sin(th)];      
     end 
     o.grid.freeGeoPts(bd,[x;y], bdsi);  % see grid2DHU  
     o.fem=lagrange12D;         
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

