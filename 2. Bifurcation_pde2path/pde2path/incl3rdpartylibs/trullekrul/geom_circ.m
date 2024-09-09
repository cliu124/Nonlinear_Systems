function dist = geom_circ(x_,r,c) %or sphere
if nargin == 1
	r=1;
end;
if nargin < 3
	c=repmat(0,1,size(x_,2));
end;
x_ = x_-repmat(c,size(x_,1),1);
rr = sqrt(sum(x_.^2,2)+eps);
dir_r = -x_./repmat(rr,1,size(x_,2)); 
dir_r(rr < eps,:) = 0;
dist = [(r-rr) ...
          dir_r];