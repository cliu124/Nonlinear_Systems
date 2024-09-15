function p=gentau(varargin)
% GENTAU: generate predictor from kernel com in c(q)swibra
% 
% p=gentau(p,v)
% p=gentau(p,v,newdir)
% p=gentau(p,v,newdir,lam')
%
% v=coefficients of tangents stored in p.mat.ker, 
% pred normalized wrt xi-norm, stored in p.tau, plotted in figure 6

lamdot=0; fname=[]; 
if isstruct(varargin{1}); p=varargin{1}; v=varargin{2}; 
 if nargin>2; fname=varargin{3}; end; 
 if nargin>3; lamdot=varargin{4}; end ; 
end
if ~isempty(fname); p=setfn(p,fname); end 
p.file.count=0; p.fuha.savefu(p); p.file.count=1;
t=0; ker=p.mat.ker; 
for i=1:size(v(:),1); t=t+ker(:,i)*v(i); end
t=t/xinorm(t,p.sol.xi,0,0); p.tau=t;
plotsolu(p,[p.tau; p.u(p.nu+1:end)],6,p.plot.pcmp,p.plot.pstyle); 
if lamdot~=0; p.tau(end)=lamdot; end 