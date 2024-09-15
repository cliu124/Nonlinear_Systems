function [x,p]=blss(A,b,p)
% BLSS: BorderedLinearSystemSolver
%
%  x=blss(A,b,p)
%
% See also bss, lss, bellss, bel, belpolss
x=A\b; 
