%% 
% HOPLOT: basic plotting routine for Hopf orbits
%
%  hoplot(p,wnr,cnr,varargin)
function hoplot2(p,wnr,cnr,varargin)
if nargin>3; aux=varargin{1}; 
else; try aux=p.hopf.plot; catch; aux=[]; end 
end 
skip=aux(1); tl=p.hopf.tl; ind=1:skip:tl; 
z=p.hopf.y(:,ind); tv=p.hopf.t(ind); T=p.hopf.T; 
[po,tr,ed]=getpte(p); ndim=size(po,1); 
figure(wnr); clf 
if ~isfield(p,'x0i'); p.x0i=1; end 
sol.x=T*tv; sol.y=z; xtplot(p,sol,wnr,cnr,[15 30],[]); 
set(gca,'FontSize',p.plot.fs); 
end 

