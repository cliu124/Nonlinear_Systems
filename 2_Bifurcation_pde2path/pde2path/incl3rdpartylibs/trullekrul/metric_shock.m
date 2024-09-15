function metric = metric_shock(x_, width, angle, pos, dxy, dxymin);
onesN = ones(size(x_,1),1); zerosN = zeros(size(x_,1),1);
x_ = x_-repmat(pos,size(x_,1),1);
eigL = @(x_) [1./(1+(dxymin/dxy-1)*exp(-(cos(angle)*x_(:,1)+sin(angle)*x_(:,2)).^2/width^2)) repmat(onesN,1,size(x_,2)-1)]/dxy;
v1xn = @(x_) cos(angle).*onesN;
v1yn = @(x_) sin(angle).*onesN;
if size(x_,2) == 3
v1zn = @(x_) zerosN;
vek = [v1xn(x_) v1yn(x_) v1zn(x_)]; %vek = vek./repmat(sqrt(sum(vek.^2,2)),1,3);
eigR = metric_vec2rthbss(vek);
else
eigR = [v1xn(x_) v1yn(x_) -v1yn(x_) v1xn(x_)];
end;
metric = analyt_prod(analyt_fulleig(eigL(x_)),eigR);