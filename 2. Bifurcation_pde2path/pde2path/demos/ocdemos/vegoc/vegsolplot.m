function vegsolplot(p,v,pfak,fn) % just to make the script shorter
s1=p.oc.s1; sol=p.cp; zdia=vegdiagn(p,15,pfak,fn); 
load('vegcm.asc'); load('watcm.asc');
% colormaps, saved via  save('watcm.asc','watcm','-ascii'); etc 
psol3D(s1,sol,1,1,[]); view(v); zlabel('v'); ylabel('t'); xlabel('x'); colormap(vegcm); % plot v 
psol3D(s1,sol,2,2,[]); view(v); zlabel('w'); ylabel('t'); xlabel('x'); colormap(watcm); % plot w 
psol3D(s1,sol,4,5,[]); view(v); zlabel('E'); ylabel('t'); xlabel('x'); colormap hot; % plot E 
