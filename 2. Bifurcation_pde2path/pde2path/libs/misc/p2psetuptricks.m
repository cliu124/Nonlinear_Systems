%%
% p2psetuptricks: various tricks to fix unexpected matlab behavior
% 
% organized in cells for quick calls of specific (often version dependent) 
% fixes; feel free to add stuff here which works for you
%% parpool startup failure (v>2018b) 
ps = parallel.Settings;  
ps.Pool.AutoCreate = false;
ps.Pool.IdleTimeout = Inf;
%% toolbox lost, path messed up, run p2phome/setpde2path afterwards 
restoredefaultpath
rehash toolboxcache
which -all pathdef
