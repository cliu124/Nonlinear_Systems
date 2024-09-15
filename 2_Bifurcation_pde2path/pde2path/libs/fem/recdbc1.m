function [geo, bc]=recdbc1(lx,ly,ssc)
% REDBC1: geometry and BC's for rectangle with Dirichlet BC, 1 component
%
%  [geo, bc]=recdbc1(lx,ly)
% Rectangle: [-lx,lx] times [-ly,ly] 
% ssc=stiff-spring constant 
%
% See also rec, recnbc1, gnbc, polygong
geo=rec(lx,ly); bc = gnbc(1,4,ssc,0);
