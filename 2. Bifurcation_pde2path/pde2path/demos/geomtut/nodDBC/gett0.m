function [t0,par]=gett0(varargin)
% ..
if nargin==2; p=varargin{1}; t0=varargin{2}; 
else p=loadp(varargin{1}, varargin{2});  t0=varargin{3}; 
end 
par=getaux(p); R=max(sqrt(p.X(:,2).^2+p.X(:,1).^2)); 
h0=p.u(p.nu+1); a=(2*h0*R-1).^2-1; par(4)=a; 
f1=@(t) 2*h0-cos(t)-sqrt(cos(t)^2+a); t0=fzero(@(t) f1(t),t0); 
