function [geo, bc]=recnbc2(lx,ly)
% RENBC2: geometry and b.c.'s for rectangle with Neumann BC, 2 components 
%
%  [geo, bc]=recnbc2(lx,ly)
% Rectangle: [-lx,lx] times [-ly,ly] 
%
% See also rec, recnbc1, recdbc1, gnbc, polygong
 geo=rec(lx,ly);
 bc=[[  2     2     2     2]; % Neumann BC, exported from gui 
     [0     0     0     0];
     [1     1     1     1];
     [1     1     1     1];
     [1     1     1     1];
     [1     1     1     1];
     [1     1     1     1];
     [1     1     1     1];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [48    48    48    48];
    [49    49    49    49];
    [48    48    48    48];
    [48    48    48    48];
    [49    49    49    49];
    [48    48    48    48];
    [48    48    48    48]];
end
