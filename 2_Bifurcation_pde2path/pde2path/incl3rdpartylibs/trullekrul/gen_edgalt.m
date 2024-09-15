function edgalt = gen_edgalt(edg,tri,edg2tri,bndedg_)
edgs = size(edg,1);
dim = size(tri,2);
bndedg = false(edgs,1); bndedg(bndedg_) = true;
Iedg = not(bndedg); edgs = edgs-numel(bndedg_);
edgalt = zeros(size(edg));
edgnds = zeros(dim*2,edgs);
edgnds(1:dim,:) = reshape(tri(edg2tri(Iedg,1),:)',dim,edgs);
edgnds(1+dim:2*dim,:) = reshape(tri(edg2tri(Iedg,2),:)',dim,edgs);
edgnds(edgnds==repmat(edg(Iedg,1)',2*dim,1)) = 0;
edgnds(edgnds==repmat(edg(Iedg,2)',2*dim,1)) = 0;
if dim==4
edgnds(edgnds==repmat(edg(Iedg,3)',2*dim,1)) = 0;
end;
edgnds = reshape(edgnds(edgnds~=0),2,edgs)'; edgalt = zeros(numel(Iedg),2); 
edgalt(Iedg,:) = edgnds;