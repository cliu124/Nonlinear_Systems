function dist = geom_cyl(x_,r,poscr) %we do not take ends into account, i.e. only use with subtraction
if nargin~=2
	x_ = x_ - repmat(poscr,size(x_,1),1);
end;
dist = [geom_circ(x_(:,[1 2]),r) zeros(size(x_,1),1)];