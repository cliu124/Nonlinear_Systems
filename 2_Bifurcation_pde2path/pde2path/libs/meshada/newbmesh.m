function p=newbmesh(p)
% NEWBMESH: inquires to build new base mesh 
%
%  p=newbmesh(p)
%
% See also stanmesh, newmesh
q=p; 
msel=0; msel=asknu('Delaunay (0) or point-mesh(1)?',msel);
if(msel==0)
    xmin=min(p.mesh.p(1,:)); xmax=max(p.mesh.p(2,:)); 
    ymin=min(p.mesh.p(2,:)); ymax=max(p.mesh.p(2,:)); 
    hmax=min((xmax-xmin)/20,(ymax-ymin)/20);
    hmax=asknu('hmax:', hmax); q=stanmesh(q,hmax);
else nx=20;ny=20; nx=asknu('nx', nx);ny=asknu('ny', ny);
    q=stanmesh(q,nx,ny);
end
p.mesh.bp=q.points; p.mesh.be=q.edges; p.mesh.bt=q.tria; p.mesh.maxt=2*q.nt; 

