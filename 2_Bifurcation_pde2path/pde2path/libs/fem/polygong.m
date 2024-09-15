function g = polygong(varargin)
% Polygong: returns a decomposed geometry matrix for a polygon object
%
%  g = polygong(varargin)
% g = polygong(pointlist) or g = polygong(xlist,ylist) by Pruefert
% Arguments can be
% (i) a pointlist of the form [x1 y1 x2 y2 ...]
% (ii) the koordinates of y and y component [x1 x2 ...] [y1 y2 ...]
% Examples
% (i) g =  polygong([0 1 0],[0 0 1]) triagle
% (ii) g =  polygong([0 0 1 0 0 1]) the same as (i)
% g = polygong(cos(linspace(0,2*pi,80)),sin(linspace(0,2*pi,80)))
% approximation on a circle by a "octadecagon" ;-) it is NOT the same as
% circleg
%
% NOTE the points of the polygon should be aranged in  counter clockwise
% order
%
% See also gnbc, rec, recdbc1, recnbc1, recnbc2
switch length(varargin)
    case 1
        pointlist = varargin{1};
        if mod(length(pointlist),2) == 1,
            error('number of x and y coordinates must be the same')
        end
    case 2
        if ~(length(varargin{1})==length(varargin{2})),
            error('xlist and ylist must have the same length')
        end
        nop = length(varargin{1});
        pointlist(1:2:2*nop) = varargin{1};
        pointlist(2:2:2*nop) = varargin{2};
    otherwise
        error('wrong number of input arguments')
end


% remove endpoint if it is also starting point
if norm(pointlist(1:2)-pointlist(end-1:end))<1e-5,
    pointlist=pointlist(1:end-2);
end
nop = length(pointlist)/2;
nokoords = length(pointlist);
xpoints = pointlist(1:2:nokoords);
ypoints = pointlist(2:2:nokoords);
g = 2*ones(1,nop);
g = [g;xpoints;xpoints(2:end),xpoints(1)];
g = [g;ypoints;ypoints(2:end),ypoints(1)];
g = [g;ones(1,nop);zeros(1,nop)];
