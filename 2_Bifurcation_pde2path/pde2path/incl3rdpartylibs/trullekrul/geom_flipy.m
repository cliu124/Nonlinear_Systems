function dist = geom_flipy(x_,func,yval)
dist = func([x_(:,1) repmat(2.*yval,size(x_,1),1)-x_(:,2)]);
dist(:,3) = -dist(:,3);