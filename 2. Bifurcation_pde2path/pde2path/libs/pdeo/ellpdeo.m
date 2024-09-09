classdef ellpdeo< pde  
% ellpdeo: ellipse pdeo
%  pdeo=ellpdeo(R,nr,nphi,ell)
%
%  R=radius, nr,nphi=#of discret.points, ell=ellipticity 

methods(Access = public)
  function o=ellpdeo(R,nr,nphi,ell) % constructor 
     o.grid=grid2D; r=linspace(1/nr,1,nr); r=r.^1.1; r=R*r; 
     x=[]; y=[]; 
     for i=1:nr; 
       th=linspace(0,2*pi,i*nphi); th=th(1:end-1); 
       x=[x ell*r(i)*cos(th)]; y=[y r(i)*sin(th)]; 
     end
     o.grid.pts2dom([x;y]); % use point2domain meshing 
     o.fem=lagrange12D; 
  end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

