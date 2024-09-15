function obj=setidssq(obj)
% setidssq: set the IDs (boundary-segment number) for a 2D rectangle
xmin=min(obj.p(1,:)); xmax=max(obj.p(1,:)); 
ymin=min(obj.p(2,:)); ymax=max(obj.p(2,:)); 
%xmin(ones(1,3)), xmax(ones(1,3)), ymax(ones(1,3)), 
tol=1e-4; 
for k=1:obj.nEdges
    A=obj.p(:,obj.e(1:2,k));% A, xmin(ones(1,3)), pause 
    if max(abs(A(1,:)-xmin(ones(1,2))))<tol; obj.e(5,k)=4;
    elseif max(abs(A(1,:)-xmax(ones(1,2))))<tol; obj.e(5,k)=2;
    elseif max(abs(A(2,:)-ymin(ones(1,2))))<tol; obj.e(5,k)=1;
  elseif max(abs(A(2,:)-ymax(ones(1,2))))<tol; obj.e(5,k)=3;
    else obj.e(5,k)=1;
  end
  obj.e(6,k)=1; obj.e(7,k)=0;
end      