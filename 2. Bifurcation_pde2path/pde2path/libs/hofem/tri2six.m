function q=tri2six(q); %(tri,po,q)
% tri2six: convert (OOPDE) 3-node triangulation to FSELIB 6-node-triangles 
% 
% store new triangulation in q.tri; boundary-info not yet treated
% if q.t2sisw=1, then also interpolate u and tau to new mesh 
tri=q.pdeo.grid.t'; po=q.pdeo.grid.p'; 
nt=size(tri,1); xp=po(:,1); yp=po(:,2);  
try; npo=q.np; nuo=q.nu; catch; npo=1; nuo=1; end 
for i=1:nt % add all edge midpoints to triangulation 
   x(i,1)=xp(tri(i,1));  y(i,1)=yp(tri(i,1)); 
   x(i,2)=xp(tri(i,2));  y(i,2)=yp(tri(i,2));  
   x(i,3)=xp(tri(i,3));  y(i,3)=yp(tri(i,3)); 
   x(i,4)=0.5*(x(i,1)+x(i,2)); y(i,4)=0.5*(y(i,1)+y(i,2)); 
   x(i,5)=0.5*(x(i,2)+x(i,3)); y(i,5)=0.5*(y(i,2)+y(i,3)); 
   x(i,6)=0.5*(x(i,3)+x(i,1)); y(i,6)=0.5*(y(i,3)+y(i,1));
end 
% six nodes of the first element are entered manually
p(1,1)=x(1,1); p(1,2)=y(1,1); p(2,1)=x(1,2); p(2,2)=y(1,2); p(3,1)=x(1,3); p(3,2)=y(1,3); 
p(4,1)=x(1,4); p(4,2)=y(1,4); p(5,1)=x(1,5); p(5,2)=y(1,5); p(6,1)=x(1,6); p(6,2)=y(1,6); 
c(1,1)=1; c(1,2)=2; c(1,3)=3; c(1,4)=4; c(1,5)=5; c(1,6)=6;  
ng=6;
% loop over further elements, Iflag=0 will signal a new global node
eps=0.000001;  % taken from FSElib 
for i=2:nt        % loop over elements
 for j=1:6          % loop over element nodes
 Iflag=0;
 for k=1:ng
  if(abs(x(i,j)-p(k,1)) < eps)
   if(abs(y(i,j)-p(k,2)) < eps)
     Iflag=1;    % the node has been recorded previously
     c(i,j)=k;   % the jth local node of element i is the kth global node
   end
  end
 end
 if(Iflag==0)  % record the node
   ng=ng+1; p(ng,1)=x(i,j); p(ng,2)=y(i,j);
   c(i,j)=ng;   % the jth local node of element is the new global node
 end
 end
end
pold=q.pdeo.grid.p'; q.pdeo.grid.pts2dom(p'); % convert new mesh to OOPDE-mesh 
q.np=size(p,1); q.nu=q.nc.neq*q.np; q.nt=q.pdeo.grid.nElements; 
q.hofem.tri=c; % store 6-node-mesh 
if q.hofem.t2sinterpol==1; % interpol soln to new mesh 
   ui=zeros(q.nc.neq*q.np,1); taui=ui; xo=pold(:,1); yo=pold(:,2); xn=p(:,1); yn=p(:,2);   
   par=q.u(nuo+1:end); np=q.np; lamd=q.tau(nuo+1:end);  %size(xo), size(xn), np,npo
   for i=1:q.nc.neq
     ui((i-1)*np+1:i*np)=p2interpol(xn,yn,q.u((i-1)*npo+1:i*npo),xo,yo,q); 
     taui((i-1)*np+1:i*np)=p2interpol(xn,yn,q.tau((i-1)*npo+1:i*npo),xo,yo,q); 
   end
   q.u=[ui;par]; q.tau=[taui;lamd]; 
end
