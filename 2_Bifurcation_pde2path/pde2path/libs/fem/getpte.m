function [po,tr,e]=getpte(p)
% GETPTE: get points, triangles (elements), edges from p
%  [po,tr,e]=getpte(p)
try po=p.pdeo.grid.p; tr=p.pdeo.grid.t; e=p.pdeo.grid.e; 
catch; po=p.mesh.p; tr=p.mesh.t; e=p.mesh.e; 
end
