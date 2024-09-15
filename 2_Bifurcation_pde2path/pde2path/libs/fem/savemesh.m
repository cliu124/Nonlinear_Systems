function p=savemesh(p)
% SAVEMESH: save mesh (needed if p.file.msave=0, i.e., mesh not saved w. points)
%
% see also meshada, meshref
if isdir(p.file.mdir)==0; mkdir(p.file.mdir); end
m=p.mesh; ptname=strrep(p.file.pname, '/','_'); 
mname=[p.file.mdir '/' ptname mat2str(p.file.count) 'mesh']; 
save(mname,'m'); p.file.mname=mname;
