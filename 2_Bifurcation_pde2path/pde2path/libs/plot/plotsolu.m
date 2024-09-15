function plotsolu(p,u,wnr,cnr,pstyle,varargin)
% PLOTSOLU: plot component of u with plotsol
%
%  plotsolu(p,u,wnr,cnr,pstyle)
%
% cnr=component number, wnr=window number, pstyle=plot style
%
% See also plotsol
p.up=u; p.u=u; plotsol(p,wnr,cnr,pstyle,varargin{:}); 