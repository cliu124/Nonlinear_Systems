function bndedg2edg = bks_bndedg2edg(edg,nd2edg,bndedg)
if size(edg,2) == 2
	%nd2edg_ = nd2edg(bndedg(:,1),:); nd2edg_(nd2edg_==0) = 		size(edg,1)+1; edg_ = [edg(:,2); 0];
	%%I = edg_(nd2edg_) == repmat(bndedg(:,2),1,size(nd2edg,2));
	%I = reshape(edg_(nd2edg_),size(nd2edg_)) == repmat(bndedg(:,2),1,size(nd2edg,2));
	%nd2edg_ = nd2edg_'; bndedg2edg = nd2edg_(I');
	
	nd2edg_ = nd2edg(bndedg(:,1),:); I = nd2edg_ ~= 0;
	edg2 = zeros(size(nd2edg_)); edg2(I) = edg(nd2edg_(I),2);
	II = edg2 == repmat(bndedg(:,2),1,size(nd2edg,2));
	nd2edg_ = nd2edg_'; bndedg2edg = nd2edg_(II');
else
	%nd2fac_ = nd2edg(bndedg(:,1),:); nd2fac_(nd2fac_==0) = size(edg,1)+1; fac2_ = [edg(:,2); 0]; fac3_ = [edg(:,3); 0];
	%%I = and(fac2_(nd2fac_) == repmat(bndedg(:,2),1,size(nd2edg,2)),fac3_(nd2fac_) == repmat(bndedg(:,3),1,size(nd2edg,2)));
	%I = and(reshape(fac2_(nd2fac_),size(nd2fac_)) == repmat(bndedg(:,2),1,size(nd2edg,2)),reshape(fac3_(nd2fac_),size(nd2fac_)) == repmat(bndedg(:,3),1,size(nd2edg,2)));
	%nd2fac_ = nd2fac_'; bndedg2edg = nd2fac_(I');
	
	nd2fac_ = nd2edg(bndedg(:,1),:); I = nd2fac_ ~= 0;
	fac2 = zeros(size(nd2fac_)); fac2(I) = edg(nd2fac_(I),2);
	fac3 = zeros(size(nd2fac_)); fac3(I) = edg(nd2fac_(I),3);
	II = and(fac2 == repmat(bndedg(:,2),1,size(nd2edg,2)),fac3 == repmat(bndedg(:,3),1,size(nd2edg,2)));
	nd2fac_ = nd2fac_'; bndedg2edg = nd2fac_(II');
end;