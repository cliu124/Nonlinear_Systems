function obj=setidsbar(obj)
% setidsbar: set the IDs (boundary-segment number) for a 3D bar
xmin=min(obj.p(1,:)); xmax=max(obj.p(1,:)); 
ymin=min(obj.p(2,:)); ymax=max(obj.p(2,:)); 
zmin=min(obj.p(3,:)); zmax=max(obj.p(3,:));
%xmin(ones(1,3)), xmax(ones(1,3)), ymax(ones(1,3)), 
tol=1e-4; %xmin, xmax,ymin, ymax,zmin, zmax
obj.e(5,:)=0; 
for k=1:obj.nEdges
    A=obj.p(:,obj.e(1:3,k));% A, xmin(ones(1,3)), pause 
    if max(abs(A(1,:)-xmin(ones(1,3))))<tol; obj.e(5,k)=2;
    elseif max(abs(A(1,:)-xmax(ones(1,3))))<tol; obj.e(5,k)=4;
    elseif max(abs(A(2,:)-ymin(ones(1,3))))<tol; obj.e(5,k)=3;
    elseif max(abs(A(2,:)-ymax(ones(1,3))))<tol; obj.e(5,k)=5;
    elseif max(abs(A(3,:)-zmin(ones(1,3))))<tol; obj.e(5,k)=1;
    elseif max(abs(A(3,:)-zmax(ones(1,3))))<tol; obj.e(5,k)=6;                
    end
end      
ids=obj.e(5,:); isz=ids==0; 
if any(isz); fprintf('some ids=0 in setidsbar, setting to 1\n'); 
  obj.e(5,isz(:))=1; 
  end 

