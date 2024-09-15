function p=stanmesh(p,varargin)
% STANMESH: generate mesh, pdetoolbox setting
%
%  p=stanmesh(p,varargin)
%  1 aux arg. hmax:  max triangle side lenght hmax 
%  2 aux arg. nx,ny: regular nx x ny mesh, only for rectangles
%
%  See also setmesh, jigglemesh, setbmesh
if(nargin<3) 
  hmax=varargin{1};[p.mesh.p,p.mesh.e,p.mesh.t]=initmesh(p.mesh.geo,'Hmax',hmax);
  p.mesh.p=jigglemesh(p.mesh.p,p.mesh.e,p.mesh.t,'opt','off','iter',10);
else nx=varargin{1};ny=varargin{2};
    [p.mesh.p,p.mesh.e,p.mesh.t]=poimesh(p.mesh.geo,nx,ny); 
    if p.mesh.sympoi==1; 
      m=p.mesh; [p1,e1,t1]=refinemesh(m.geo,m.p,m.e,m.t,'longest') ;
      p.mesh.p=p1; p.mesh.e=e1; p.mesh.t=t1;
    end
end
figure(p.plot.ifig);pdemesh(p.mesh.p,p.mesh.e,p.mesh.t); axis tight; % plot mesh
p.np=size(p.mesh.p,2); p.nu=p.np*p.nc.neq; p.mesh.nt=size(p.mesh.t,2); 
p.mesh.maxt=2*p.mesh.nt; 
