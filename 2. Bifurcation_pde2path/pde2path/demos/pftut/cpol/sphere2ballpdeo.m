%% 
% augment points from sphere by points in ball and call Delauney 
% -> first p0.nu points in p correspond to points in p0
% bulk and surf coupling at these points! 
classdef sphere2ballpdeo< pde  
methods(Static, Access=public)
  function o=sphere2ballpdeo(p0,R,h) % constructor 
    gr=grid3D; 
    pde0=p0.pdeo; gr0=pde0.grid;  x0f=gr0.p(1,:)'; y0f=gr0.p(2,:)'; 
    x0=p0.mat.drop*x0f; y0=p0.mat.drop*y0f; x0=x0'; y0=y0'; % drop periodic points 
    x=[cos(y0).*cos(x0) 0 0]; % map to unit sphere and add N and S(outhpole) 
    y=[cos(y0).*sin(x0) 0 0]; 
    z=[sin(y0) 1 -1]; 
    p=R*[x; y; z]';         % first p0.nu points on sphere of radius R 
    pfix=p; box=[-1 1;-1 1;-1 1]; box=box'; imax=2; 
    [p,e,t]=hudistmesh(@fd,@fh, h, box, imax, pfix); 
    t=[t; ones(1,size(t,2))]; % pad tetra list with ones 
    se=size(e,2); ed=[e; zeros(2,se)]; % pad edgelist with 0 
    ed(5,:)=ones(1,se); % and assign each face with segment nr 1
    gr.p=p; gr.t=t; gr.e=ed; o.grid=gr; o.fem=lagrange13D; 
    function d=fd(p)         
        d=sqrt(sum(p.^2, 2))-(1-h/2); 
    end
    function h=fh(p, varargin)   
        np=size(p,1); h=ones(np,1);
    end  
  end
  
end
methods(Access=protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end

