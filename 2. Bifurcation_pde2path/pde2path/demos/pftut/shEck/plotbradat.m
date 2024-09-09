function plotbradat(p,wnr,xc,yc,varargin)
% PLOTBRADAT: plot branch component vs. branch component
% xc=x-axis component; yc=y-axis component
if nargin>4; aux=varargin{1}; end 
try ps=aux.ps; catch ps='-k'; end 
figure(wnr); plot(p.branch(6+xc,:),p.branch(6+yc,:),ps);
