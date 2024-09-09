function dist = geom_point(x_,p)
dr = x_-repmat(p,size(x_,1),1);
r = sqrt(sum(dr.^2,2));
I = r > eps; dr(I,:) = dr(I,:)./repmat(r(I),1,size(x_,2));
dist = [r dr];