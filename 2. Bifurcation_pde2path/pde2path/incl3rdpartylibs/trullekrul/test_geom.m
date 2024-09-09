clc; clear all; close all;
testN = 1;

switch testN
case 1
N = 31; x=linspace(-0.5,2.5,N); y = linspace(-0.1,1.1,N); [X,Y]=meshgrid(x,y);
distC = geom_circ([X(:) Y(:)],0.25,[0 0.5]); 
distR = geom_rect_([X(:) Y(:)],[2 1],[0 0]); 
%dist = distC
%dist = distR
%dist = geom_union(distR,distC);
%dist = geom_intersect(distR,distC);
dist = geom_diff(distR,distC);
dist = geom_diff(distC,distR);
Z = reshape(distR(:,1),size(X)); Zx = reshape(dist(:,2),size(X)); Zy = reshape(dist(:,3),size(X));
contourf(X,Y,Z,100); colorbar(); hold on;% quiver(X,Y,Zx,Zy,'color',[1 1 1]*0.99); 
%contour(X,Y,Z,[eps 0],'k'); 
I = abs(Z(:)) < 0.02; plot(X(I),Y(I),'.r'); hold off;
hold off; axis equal;
case 2
cnst = [1.5 0.5 0.5]; inds = [1 2 3; 1 3 2; 3 2 1];
N = 31; x=linspace(-0.5,2.5,N); y = linspace(-0.1,1.1,N); [X,Y]=meshgrid(x,y);
for i=1:3
coords = [X(:) Y(:) repmat(cnst(i),numel(X),1)];
%dist = geom_circ(coords(:,inds(i,:)),0.25,[0.5 0.5 0.5]); 
dist = geom_cyl(coords(:,inds(i,:)),[0.25, 0.25, 0.25],[0.5 0.5 0.5]); 
%distR = geom_rect([X(:) Y(:) repmat(cnst(i),numel(X),1)],[2 1 1],[0 0 0]); 
Z = reshape(dist(:,1),size(X)); 
Zx = reshape(dist(:,1+inds(i,1)),size(X)); 
Zy = reshape(dist(:,1+inds(i,2)),size(X));
%figure(i); 
%contourf(X,Y,Z,100);  hold on; %colorbar(); quiver(X,Y,Zx,Zy,'color',[1 1 1]*0.99); 
%contour(X,Y,Z,[eps 0],'k');  title(sprintf('x%0.0f=%2.2f',inds(i,3),cnst(i)))
hold off; axis equal;
end;

end;
