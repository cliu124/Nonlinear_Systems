function [x,p]=amgc(A,b,p)
% amgc: amg by chen
%
%  [x,p]=lss(A,b,p)
%
% See also blss, ilss
x=amg(A,b,p.opt); 
