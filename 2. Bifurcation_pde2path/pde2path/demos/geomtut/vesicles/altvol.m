function V=altvol(p) % alternate(?) volume 
VF=permute(reshape(p.X(p.tri,:),[size(p.tri) 3]),[3 1 2]);
V=1/6*abs(sum(dot(cross(VF(:,:,1),VF(:,:,2),1),VF(:,:,3),1))); 
