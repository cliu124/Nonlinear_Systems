function plotsolusf(p,u,wnr,cnr,pstyle,m,n,pos)
%PLOTSOLU: plot component of u with plotsol 
%   plotsolu(p,u,wnr,cnr,pstyle)
% cnr=component number, wnr=window number, pstyle=plot style
%
% See also plotsol
p.u=u; plotsolsf(p,wnr,cnr,pstyle,m,n,pos); 