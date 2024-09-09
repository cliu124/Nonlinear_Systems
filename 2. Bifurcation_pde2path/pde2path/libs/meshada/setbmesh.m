function p=setbmesh(p)
% SETBMESH: (re)set basemesh for meshadac to current mesh.
%
%  p=setbmesh(p) 
%
% See also meshadac, stanparam.
try; p.mesh.bp=p.mesh.p; p.mesh.be=p.mesh.e; p.mesh.bt=p.mesh.t; p.mesh.maxt=2*p.mesh.nt; 
catch; p.mesh.bp=p.pdeo.grid.p; p.mesh.be=p.pdeo.grid.e; p.mesh.bt=p.pdeo.grid.t; 
end
