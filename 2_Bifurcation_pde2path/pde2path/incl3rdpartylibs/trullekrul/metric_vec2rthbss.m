function rthbss = metric_vec2rthbss(vek)
N = size(vek,1);
[tmp,Cd] = min(abs(vek),[],2);
v2 = zeros(size(vek)); v2((1:N)'+(Cd-1)*N) = 1;
v3 = cross(vek,v2,2);
v3 = v3./repmat(sqrt(sum(v3.^2,2)),1,3);
v2 = cross(vek,v3,2);
rthbss = [vek v2 v3];