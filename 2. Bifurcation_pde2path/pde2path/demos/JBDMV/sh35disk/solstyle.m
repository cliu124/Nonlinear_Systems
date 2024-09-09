function solstyle(txt,varargin) % modify soln plots 
nolti; colormap parula; title(''); colorbar off; h=colorbar('southoutside'); 
%h.Ticks=[-0.5 0.5]; % for low ampl.soln plots, adapt ticks of colorbar by hand
if nargin==1; text(9,12,txt,'fontsize',16);
else text(9,12,txt,'fontsize',16,'color',varargin{1});
end 