p=loadp('ws','pt0'); plotsol(p,1,1,3); myticks(['org, n_p=' mat2str(p.np)]); view(v); 
%plotsol(p,6,1,4); noticks(['org, n_p=' mat2str(p.np)]); zticks([-2 0]); view(v); 
%%
p=loadp('wsr','pt0'); plotsol(p,1,1,3); myticks(['refinement, n_p=' mat2str(p.np)]); view(v); 
%%
p=loadp('wsrc','pt0'); plotsol(p,1,1,3); myticks(['extra coarsening, n_p=' mat2str(p.np)]); view(v); 