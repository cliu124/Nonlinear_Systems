classdef spherepdeo < pde  
% spherepdeo: (preimage-)domain and mesh for a sphere 
% modification of stanpdeo2D to have coarser grids near poles  
%
%  pdeo=spherepdeo(lx,ly,n1,n2,ref,yexp)
%
% lx=dom.size (preimage) in x, usually lx=pi; 
% ly=dom.size (preimage) in y, usually lx=pi/2-del to exclude poles  
% n1,n2=#points in x and y
% ref=#initial refinement steps, 
% yexp=scaling in elevation (make mesh coarser near poles
% 
methods(Access = public)
  function o=spherepdeo(lx,ly,n1,n2,ref,yexp)     
     t1=linspace(-1,1,n2); t1=sign(t1).*abs(t1).^yexp; 
     x=[]; y=[]; P=[]; 
     for i=1:n2
         nx=n2+n1-round(abs(n2*t1(i))); 
         xa=linspace(-lx,lx,nx); ya=ly*t1(i)*ones(1,nx); 
         x=[x xa]; y=[y ya]; P=[P [xa; ya]]; 
     end
     P=P'; P=unique(P,'rows','stable'); P=P'; 
     dt=delaunayTriangulation(P'); po=dt.Points'; e=dt.freeBoundary'; t=dt.ConnectivityList';        
     t=[t; ones(1,size(t,2))]; 
     se=size(e,2); ed=[e; zeros(3,se)]; ed(5,:)=ones(1,se);
     ed(4,:)=ones(1,se); ed(3,:)=ones(1,se); 
     grid=grid2DHU; grid.p=po; grid.t=t; grid.e=ed; 
     o.grid=grid; o.fem=lagrange12D;
     for i=1:ref; o.grid.refineMesh; end;
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
