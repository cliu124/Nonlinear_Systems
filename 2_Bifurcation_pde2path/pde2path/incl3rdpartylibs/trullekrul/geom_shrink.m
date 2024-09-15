function dist = geom_shrink(x_,func,shrink,move)
if nargin == 3
move = zeros(1,size(x_,2));
end;
dist = func(x_*shrink+repmat(move,size(x_,1),1)*(1.-shrink));
dist(:,1) = dist(:,1)/shrink;