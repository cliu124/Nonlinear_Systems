function [x,p]=lss(A,b,p)
% LSS: default LinearSystemSolver, i.e. \ 
%
%  [x,p]=lss(A,b,p)
%
% See also blss, ilss
x=A\b; 
