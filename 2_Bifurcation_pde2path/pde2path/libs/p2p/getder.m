function [Gu,Glam]=getder(p,u,r)
% GETDER: return jacobian Gu and derivative Glam depending on p.sw.jac 
%
%  [Gu,Glam]=getder(p,u,r)
%
% See also getGu, getGlam, stanparam
Gu=getGu(p,u,r);Glam=getGlam(p,u,r); 
