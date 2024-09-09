classdef secpdeo1< pde  
% secpdeo1: sector pdeo contructor, simple meshing 
%
% pdeo=secpdeo(r,phi,nr,nphi) 
%
% r=radius, phi=half-angle, nr,nphi=#point in r and in phi 
methods(Access = public)
  function o=secpdeo1(r,phi,nr,nphi) % constructor 
     o.grid=grid2D;  s=linspace(-phi,phi,nphi); s=s(2:end-1); 
     bd1=linspace(0,r,nr); bda=bd1(1:end);  bdb=r-bd1(1:end-1);      
     bd=[bda*cos(-phi) r*cos(s) cos(phi)*bdb; 
         bda*sin(-phi) r*sin(s) sin(phi)*bdb]; 
     o.grid.freeGeometry(bd); o.fem=lagrange12D;         
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

