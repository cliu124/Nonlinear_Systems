function h = huafixed(varargin)
% huafixed: hack of MATLAB's annotation, use in plotbra with fancy=2
% h=huafixed(varargin)
%ANNOTATION_FIXED works exactly as the standard MATLAB function ANNOTATION,
%   except the position/width/heigth of the annotation object are given in
%   coordinates of the current axes. The annotation object is pinned to the
%   current axes so it's position doesn't change on resize.
switch varargin{1};
    case {'arrow','textarrow','doublearrow','line'}
    if verLessThan('matlab','8.4.0') % execute code for R2014a or earlier
        warning('Old Matlab. Style fancy 2 might not work. Consider to change to fancy 0,1.');
        x = varargin{2}; y = varargin{3}; % x,y 
        si=1; %0*max(x(2)-x(1), y(2)-y(1)); 
        pos = [x(1) y(1) si*(x(2)-x(1)) si*(y(2)-y(1))];
        h = annotation(varargin{1}, varargin{4:end});
        set(h,'parent',gca); set(h,'position',pos);
    else
    % execute code for R2014b or later
        x = varargin{2}(1:2); y = varargin{3}(1:2);
        xmi=varargin{2}(3); ymi=varargin{3}(3);
        xext=varargin{2}(4); yext=varargin{3}(4);
        xrel=(x-xmi)./xext; yrel=(y-ymi)./yext;
        h = annotation(varargin{1},xrel,yrel,varargin{4:end});
        h.pinAtAffordance(1);  h.pinAtAffordance(2); % not in 2013, HU
    end
    case {'rectangle','ellipse','textbox'} % probably won't work
        pos = varargin{2};
        h = annotation(varargin{1},varargin{3:end});
        set(h,'parent',gca); % h transforms with axis 
        set(h,'position',pos);
end
end 
