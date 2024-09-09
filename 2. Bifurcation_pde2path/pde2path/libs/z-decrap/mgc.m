function [x,p]=mgc(A,b,p)
% amgc: amg by chen
%
%  [x,p]=lss(A,b,p)
%
% See also blss, ilss
nd=size(p.pdeo.grid.p,1);% nd, pause 
x=mg(A,b,p.pdeo.grid.t(1:nd+1,:)',p.opt); 
