function screenlayout(varargin)
% SCREENLAYOUT: open, clear and position windows for sol, branch and info 
%
%  screenlayout(p)
%
% See also stanparam
set(0,'Units','pixels'); scnsize = get(0,'ScreenSize'); scnw=scnsize(3); scnh=scnsize(4);
if nargin==1; p=varargin{1}; 
figure(p.plot.pfig); clf(p.plot.pfig); set(figure(p.plot.pfig),'Position', [scnw-650 scnh-400 300 300]); %  profile figure
figure(p.plot.brfig); clf(p.plot.brfig); set(figure(p.plot.brfig),'Position', [scnw-300 scnh-400 300 300]); %  branch figure
figure(p.plot.ifig); clf(p.plot.ifig); set(figure(p.plot.ifig),'Position', [scnw-300 scnh-800 300 300]); %  info figure
else 
set(figure(1),'Position', [scnw-650 scnh-400 300 300]); %  profile figure
set(figure(2),'Position', [scnw-300 scnh-400 300 300]); %  branch figure
set(figure(6),'Position', [scnw-300 scnh-800 300 300]); %  info figure
end 