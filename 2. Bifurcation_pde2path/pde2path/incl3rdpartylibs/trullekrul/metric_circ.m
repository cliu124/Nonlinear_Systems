function metric = metric_circ(x_, width, radius, pos, dxy, dxymin);
onesN = ones(size(x_,1),1); zerosN = zeros(size(x_,1),1);
x_(:,[1 2]) =  x_(:,[1 2]) - repmat(pos,size(x_,1),1);
eigL = @(x_) [1./(1+(dxymin/dxy-1)*exp(-(sqrt(sum(x_(:,[1 2]).^2,2))-radius).^2/width^2)) repmat(onesN,1,size(x_,2)-1)]/dxy;
v1xn = @(x_) x_(:,1)./sqrt(sum(x_(:,[1 2]).^2,2));
v1yn = @(x_) x_(:,2)./sqrt(sum(x_(:,[1 2]).^2,2));
if size(x_,2) == 3
	v1zn = @(x_) zerosN;
	vek = [v1xn(x_) v1yn(x_) v1zn(x_)]; %vek = vek./repmat(sqrt(sum(vek.^2,2)),1,3);
	eigR = metric_vec2rthbss(vek);
else
	eigR = [v1xn(x_) v1yn(x_) -v1yn(x_) v1xn(x_)];
end;

metric = analyt_prod(analyt_fulleig(eigL(x_)),eigR); 