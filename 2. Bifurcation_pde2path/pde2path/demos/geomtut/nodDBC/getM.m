function M=getM(p,varargin) % getM: mass matrix M for p.X, or for X=varargin 
try Xcont=p.sw.Xcont; catch Xcont=0; end 
if Xcont
if nargin==1; X=p.X; else X=varargin{1}; end; 
M=massmatrix(X,p.tri,'voronoi');
%M=massmatrix(X,p.tri,'full');
else M=p.mat.M; 
end 