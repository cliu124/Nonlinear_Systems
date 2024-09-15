function vegsolplot(sol,v,pfak,fn) % just to make the script shorter
global s1; 
zdia=vegdiagn(sol,15,pfak,fn); 
load('vegcm.asc'); load('watcm.asc');
% colormaps, saved via  save('watcm.asc','watcm','-ascii'); etc 
psol3Db(s1,sol,1,1,[]); view(v); zlabel('v'); ylabel('t'); colormap(vegcm); % plot v 
psol3Db(s1,sol,2,2,[]); view(v); zlabel('w'); ylabel('t'); colormap(watcm); % plot w 
psol3Db(s1,sol,4,5,[]); view(v); zlabel('E'); ylabel('t'); colormap hot; % plot E 
