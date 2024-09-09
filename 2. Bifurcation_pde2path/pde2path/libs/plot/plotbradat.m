function h=plotbradat(p,wnr,xc,yc,varargin)
% PLOTBRADAT: plot branch component vs. branch component
%
%  plotbradat(p,wnr,xc,yc)
% xc=x-axis component
% yc=y-axis component
%
% See also plotbra, bradat, stanbra
figure(wnr);
if nargin==4;  h=plot(p.branch(xc,:),p.branch(yc,:),'*-');
else  h=plot(p.branch(xc,:),p.branch(yc,:),varargin{1});
end 
