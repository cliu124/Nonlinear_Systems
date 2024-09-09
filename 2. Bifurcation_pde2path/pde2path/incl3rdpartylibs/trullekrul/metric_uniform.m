function metric = metric_uniform(x_, dxy);
n = size(x_,1);
metric = [1./dxy(1) 0 1./dxy(2)];
if size(x_,2) == 3
	metric = [metric 0 0 1./dxy(3)];
end;
metric = repmat(metric,n,1);