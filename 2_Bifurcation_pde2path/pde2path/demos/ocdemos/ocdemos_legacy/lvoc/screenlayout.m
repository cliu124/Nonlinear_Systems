%%
% SCREENLAYOUT: open, clear and position profile, branch and info windows
%
%  screenlayout(p)
%
% See also stanparam
function screenlayout(p)
set(0,'Units','pixels'); scnsize = get(0,'ScreenSize'); scnw=scnsize(3); scnh=scnsize(4);
figure(p.plot.pfig); clf(p.plot.pfig); set(figure(p.plot.pfig),'Position', [scnw-750 scnh-400 400 300]); %  profile figure
figure(p.plot.brfig); clf(p.plot.brfig); set(figure(p.plot.brfig),'Position', [scnw-300 scnh-400 300 300]); %  branch figure
figure(p.plot.ifig); clf(p.plot.ifig); set(figure(p.plot.ifig),'Position', [0 scnh-400 300 300]); %  info figure
end