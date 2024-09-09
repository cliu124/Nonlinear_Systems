function varargout=ptselect(pdir)
% ptselect: read contents of dir and ask user (with default) 
%
%  pt=ptselect(dir) or ptselect(dir)
dn=[pdir '/pt*.mat']; pts=dir(dn);
fprintf(['Directory ' pdir ': available regular points are (with lambda values)\n']); 
labs=getlabs(pdir); 
[slabs,idx]=sort(labs);  
for j=1:length(slabs)
  cj=idx(j); p=loadp(pdir,pts(cj).name(1:end-4)); 
  fprintf([pts(cj).name(1:end-4) ', %4.4f     '], getlam(p)); 
  if mod(j,3)==0 fprintf('\n'); end 
end
dn=[pdir '/bpt*.mat']; pts=dir(dn); 
dn=[pdir '/hpt*.mat']; pts=[pts dir(dn)]; 
fprintf('\nOther points: ');
for j=1:length(pts)
  p=loadp(pdir,pts(j).name(1:end-4)); 
  fprintf([pts(j).name(1:end-4) ', %4.4f    '], getlam(p)); 
  if mod(j,3)==0 fprintf('\n'); end 
end
fprintf('\n'); 
if nargout>0; maxpt=char(['pt' mat2str(max(getlabs(pdir)))]);
varargout{1}=asks('choose point ', maxpt); 
end 