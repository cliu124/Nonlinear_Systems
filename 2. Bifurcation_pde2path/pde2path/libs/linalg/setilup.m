function p=setilup(p,varargin); %dtol,maxit)
% setilup:convenience function for setting AMGsolver (ilupack) parameters 
% p=setilup(p,dtol,maxit)
%
p.fuha.lss=@lssAMG; p.fuha.blss=@lssAMG; p.fuha.innerlss=@lssAMG; 
if nargin>1; dtol=varargin{1}; p.ilup.droptol=dtol; p.ilup.droptolS=dtol/10; end 
if nargin>2; p.ilup.maxit=varargin{2}; end 