function nids=rmbdtrifromlist(tri,ids) 
% rmbdtrifromlist: remove boundary triangles from list; for corsen and refine 
nr=length(ids); bdid=boundary_faces(tri); bdid=unique(bdid(:)); 
bdtri=[]; 
for i=1:nr; if ~isempty(intersect(tri(ids(i),:),bdid)); bdtri=[bdtri ids(i)]; end; end
nids=setdiff(ids,bdtri); 