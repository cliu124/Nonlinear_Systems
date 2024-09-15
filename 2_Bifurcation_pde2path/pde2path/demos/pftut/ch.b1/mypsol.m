function mypsol(dir,pt,varargin)
plotsol(dir,pt); nolti; axis 'image'; title([dir '/' pt]);
try view(varargin{1}); zticks([]); zlabel(''); catch; end; pause 