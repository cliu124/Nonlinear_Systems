function q=four2ten(q)
% four2ten: convert (OOPDE) 4-node tetras to FSELIB 10-node-tetras 
% 
% store new triangulation in q.po and q.tri; 
% if q.f2tisw=1, then also interpolate u and tau to new mesh 
tri=q.pdeo.grid.t'; po=q.pdeo.grid.p'; 
nt=size(tri,1); xp=po(:,1); yp=po(:,2); zp=po(:,3);   
try; npo=q.np; nuo=q.nu; catch; npo=1; nuo=1; end 
for i=1:nt % add all edge midpoints to triangulation 
 x(i,1)=xp(tri(i,1));  y(i,1)=yp(tri(i,1)); z(i,1)=zp(tri(i,1)); 
 x(i,2)=xp(tri(i,2));  y(i,2)=yp(tri(i,2)); z(i,2)=zp(tri(i,2)); 
 x(i,3)=xp(tri(i,3));  y(i,3)=yp(tri(i,3)); z(i,3)=zp(tri(i,3)); 
 x(i,4)=xp(tri(i,4));  y(i,4)=yp(tri(i,4)); z(i,4)=zp(tri(i,4)); 
 x(i,5)=0.5*(x(i,1)+x(i,2)); y(i,5)=0.5*(y(i,1)+y(i,2)); z(i,5)=0.5*(z(i,1)+z(i,2)); 
 x(i,6)=0.5*(x(i,2)+x(i,3)); y(i,6) =0.5*(y(i,2)+y(i,3)); z(i,6) =0.5*(z(i,2)+z(i,3)); 
 x(i,7)=0.5*(x(i,1)+x(i,3)); y(i,7)=0.5*(y(i,1)+y(i,3)); z(i,7)=0.5*(z(i,1)+z(i,3)); 
 x(i,8)=0.5*(x(i,1)+x(i,4)); y(i,8)=0.5*(y(i,1)+y(i,4)); z(i,8)=0.5*(z(i,1)+z(i,4));
 x(i,9)=0.5*(x(i,2)+x(i,4)); y(i,9)=0.5*(y(i,2)+y(i,4)); z(i,9)=0.5*(z(i,2)+z(i,4));
 x(i,10)=0.5*(x(i,3)+x(i,4)); y(i,10)=0.5*(y(i,3)+y(i,4)); z(i,10)=0.5*(z(i,3)+z(i,4));  
end 
% ten nodes of the first element are entered manually
p(1,1)=x(1,1);  p(1,2)=y(1,1);  p(1,3)=z(1,1);  
p(2,1)=x(1,2);  p(2,2)=y(1,2);  p(2,3)=z(1,2);  
p(3,1)=x(1,3);  p(3,2)=y(1,3);  p(3,3)=z(1,3);  
p(4,1)=x(1,4);  p(4,2)=y(1,4);  p(4,3)=z(1,4);  
p(5,1)=x(1,5);  p(5,2)=y(1,5);  p(5,3)=z(1,5);  
p(6,1)=x(1,6);  p(6,2)=y(1,6);  p(6,3)=z(1,6);  
p(7,1)=x(1,7);  p(7,2)=y(1,7);  p(7,3)=z(1,7);  
p(8,1)=x(1,8);  p(8,2)=y(1,8);  p(8,3)=z(1,8);  
p(9,1)=x(1,9);  p(9,2)=y(1,9);  p(9,3)=z(1,9);  
p(10,1)=x(1,10); p(10,2)=y(1,10); p(10,3)=z(1,10); 

c(1,1)=1; c(1,2)=2; c(1,3)=3; c(1,4)=4; c(1,5)=5; 
c(1,6)=6; c(1,7)=7; c(1,8)=8; c(1,9)=9; c(1,10)=10; 
ng=10;
% loop over further elements, Iflag=0 will signal a new global node
eps=0.000001;  % taken from FSElib 
for i=2:nt        % loop over elements
 for j=1:10         % loop over element nodes
 Iflag=0;
 for k=1:ng
  if(abs(x(i,j)-p(k,1)) < eps)
   if(abs(y(i,j)-p(k,2)) < eps)
    if(abs(z(i,j)-p(k,3)) < eps)  
     Iflag=1;  c(i,j)=k;   % the jth local node of element i is the kth global node
   end
   end
  end
 end
 if(Iflag==0)  % record the node
   ng=ng+1; p(ng,1)=x(i,j); p(ng,2)=y(i,j); p(ng,3)=z(i,j); 
   c(i,j)=ng;   % the jth local node of element is the new global node
 end
 end
end
% convert 10-nodes to 4-node mesh; for plotting, and, e.g., generate BC matrices
pold=q.pdeo.grid.p'; q.pdeo.grid.pts2dom(p'); q.np=size(p,1); q.nu=q.nc.neq*q.np; 
q.nt=q.pdeo.grid.nElements; q.hofem.tri=c; % store 10-node-mesh 
%return 
try f2tisw=q.hofem.t2sinterpol; catch f2tisw=0; end 
if f2tisw; % interpol soln to new mesh   
   ui=zeros(q.nc.neq*q.np,1); taui=ui; xo=pold(:,1); yo=pold(:,2); zo=pold(:,3); 
   xn=p(:,1); yn=p(:,2); zn=p(:,3);    
   par=q.u(nuo+1:end); np=q.np; lamd=q.tau(nuo+1:end);  %size(xo), size(xn), np,npo
   for i=1:q.nc.neq
     ui((i-1)*np+1:i*np)=p3interpol(xn,yn,zn,q.u((i-1)*npo+1:i*npo),xo,yo,zo,q); 
     taui((i-1)*np+1:i*np)=p3interpol(xn,yn,zn,q.tau((i-1)*npo+1:i*npo),xo,yo,zo,q); 
   end
   q.u=[ui;par]; q.tau=[taui;lamd]; 
end
