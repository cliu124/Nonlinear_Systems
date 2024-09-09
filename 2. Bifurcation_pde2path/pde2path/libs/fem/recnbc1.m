function [geo, bc]=recnbc1(lx,ly)
% RENBC1: geometry and b.c.'s for rectangle with Neumann BC, 1 component
%
%  [geo, bc]=recnbc1(lx,ly)
% Rectangle: [-lx,lx] times [-ly,ly] 
%
% See also rec, recnbc2, recdbc1, gnbc, polygong
geo=rec(lx,ly); 
bc=[[  1     1     1     1 ];
     [0     0     0     0];
     [1     1     1     1];
     [1     1     1     1];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [49    49    49    49];
    [48    48    48    48]];
end
