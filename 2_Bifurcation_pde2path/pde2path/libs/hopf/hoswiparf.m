function p=hoswiparf(dir,fname,ndir,par,ds,varargin)
% SWIPARF: switch parameter for (Hopf) solution from file (old name) 
%
%  p=hoswiparf(dir,fname,ndir,par,ds,varargin)
%
% varargin=new hopf.ilam
%
% See also swipar, loadp
p=loadp(dir,fname,ndir); p=resetc(p); p.nc.ilam=par; p.hopf.lam=getlam(p)
p.sol.ds=ds; 
if nargin>5; p.hopf.ilam=varargin{1}; end 