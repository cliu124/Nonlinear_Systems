function M=getM(p,varargin) % getM: mass matrix M for p.X, or for X=varargin 
if nargin==1; X=p.X; else X=varargin{1}; end; 
try msw=p.sw.msw; catch; msw=0; end 
if msw==0;  M=massmatrix(X,p.tri,'voronoi'); 
else  M=massmatrix(X,p.tri,'full'); 
end

