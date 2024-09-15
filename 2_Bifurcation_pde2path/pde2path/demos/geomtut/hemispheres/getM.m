function M=getM(p,varargin) % getM: mass matrix M for p.X, or for X=varargin 
if nargin==1; X=p.X; else X=varargin{1}; end; 
M=massmatrix(X,p.tri,'voronoi');