function dist = geom_rotate(x_,func,angle,origin)
if nargin == 3
origin = zeros(1,size(x_,2));
end;
x_ = x_ - repmat(origin,size(x_,1),1);
if size(x_,2) == 2
x_ = [cos(angle)*x_(:,1) + sin(angle)*x_(:,2), -sin(angle)*x_(:,1) + cos(angle)*x_(:,2)];
else
error('3D rotation not implemented');
end;
dist = func(x_+repmat(origin,size(x_,1),1));
dist(:,2:3) = [cos(angle)*dist(:,2) - sin(angle)*dist(:,3), sin(angle)*dist(:,2) + cos(angle)*dist(:,3)];