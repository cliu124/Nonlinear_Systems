function geo=rec(lx,ly)
% REC: provide rectangle-geometry in pdetoolbox-form   
%
%  geo=rec(lx,ly)
%
% Rectangle: [-lx,lx] times [-ly,ly]
%
% See also rec4
geo= [[2.0000    2.0000    2.0000    2.0000]; 
     [-lx     lx  lx  -lx]; % each colum: xstart
     [lx    lx   -lx   -lx]; %             xend
     [-ly    -ly    ly  ly]; %            ystart
     [-ly   ly   ly    -ly]; %            yend
     [1.0000    1.0000    1.0000    1.0000];
     [    0         0         0         0 ]];
end
