function M=getM(p,varargin) % getM: mass matrix M from p.X, or from aux.X 
if nargin==1; X=p.X; else X=varargin{1}; end; 
M=massmatrix(X,p.tri,'voronoi');
M=[M 0*M;0*M 0*M]; 