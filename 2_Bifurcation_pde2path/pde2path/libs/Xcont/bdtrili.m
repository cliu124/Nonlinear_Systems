function bdtri=bdtrili(tri) 
% bdtrili: used in rmdegfaces 
nr=size(tri,1); bdid=boundary_faces(tri); bdid=unique(bdid(:)); 
bdtri=[]; 
for i=1:nr; if ~isempty(intersect(tri(i,:),bdid)); bdtri=[bdtri i]; end; end