function [ne,ng,p,c,efl,gfl]=trgl6_rec(lx,ly,nx,ndiv)
% trgl6_rec: triangulate rectangle into 6-node traingles, following FSELIB 
% rec=(0,lx) x (0,ly), nx=discr.points x/discr. points y  (natural number) 
% ndiv=#refinement steps 
% somewhat obsolete since we rather convert OOPDE discretizations to
% 6-nodes 
ne=2*nx;  xp=0;   % 'ancestral discr' 
for i=0:nx; 
x(1+2*i,1)=xp+0.0; y(1+2*i,1)=0.0; efl(1+2*i,1)=0; % first element
x(1+2*i,2)=xp+1.0; y(1+2*i,2)=0.0; efl(1+2*i,2)=1;
x(1+2*i,3)=xp+1.0; y(1+2*i,3)=1.0; efl(1+2*i,3)=1;
x(1+2*i,4)=0.5*(x(1+2*i,1)+x(1+2*i,2)); y(1+2*i,4)=0.5*(y(1+2*i,1)+y(1+2*i,2)); efl(1+2*i,4)=0;
x(1+2*i,5)=0.5*(x(1+2*i,2)+x(1+2*i,3)); y(1+2*i,5)=0.5*(y(1+2*i,2)+y(1+2*i,3)); efl(1+2*i,5)=1;
x(1+2*i,6)=0.5*(x(1+2*i,3)+x(1+2*i,1)); y(1+2*i,6)=0.5*(y(1+2*i,3)+y(1+2*i,1)); efl(1+2*i,6)=0;

x(2+2*i,1)=xp+0.0; y(2+2*i,1)=0.0; efl(2+2*i,1)=0;  % second element
x(2+2*i,2)=xp+1.0; y(2+2*i,2)=1.0; efl(2+2*i,2)=1;
x(2+2*i,3)=xp+0.0; y(2+2*i,3)=1.0; efl(2+2*i,3)=1;
x(2+2*i,4)=0.5*(x(2+2*i,1)+x(2+2*i,2));y(2+2*i,4)=0.5*(y(2+2*i,1)+y(2+2*i,2));efl(2+2*i,4)=0;
x(2+2*i,5)=0.5*(x(2+2*i,2)+x(2+2*i,3));y(2+2*i,5)=0.5*(y(2+2*i,2)+y(2+2*i,3));efl(2+2*i,5)=1;
x(2+2*i,6)=0.5*(x(2+2*i,3)+x(2+2*i,1));y(2+2*i,6)=0.5*(y(2+2*i,3)+y(2+2*i,1));efl(2+2*i,6)=0;
xp=xp+1; 
end

% refinement loop
if(ndiv > 0)
for i=1:ndiv
 nm=0; % counter for new elements  (four new in each pass)
 for j=1:ne   % loop over current elements
  % assign vertex nodes to sub-elements, these will become the "new" elements
   nm=nm+1;
   xn(nm,1)=x(j,1); yn(nm,1)=y(j,1); efln(nm,1)=efl(j,1); %  first sub-element
   xn(nm,2)=x(j,4); yn(nm,2)=y(j,4); efln(nm,2)=efl(j,4);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);
   xn(nm,4)=0.5*(xn(nm,1)+xn(nm,2));yn(nm,4)=0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5)=0.5*(xn(nm,2)+xn(nm,3));yn(nm,5)=0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6)=0.5*(xn(nm,3)+xn(nm,1));yn(nm,6)=0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4)=0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4)=1; end
   efln(nm,5)=0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5)=1; end
   efln(nm,6)=0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6)=1; end

   nm=nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  second sub-element
   xn(nm,2)=x(j,2); yn(nm,2)=y(j,2); efln(nm,2)=efl(j,2);
   xn(nm,3)=x(j,5); yn(nm,3)=y(j,5); efln(nm,3)=efl(j,5);
   xn(nm,4)=0.5*(xn(nm,1)+xn(nm,2));yn(nm,4)=0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5)=0.5*(xn(nm,2)+xn(nm,3));yn(nm,5)=0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6)=0.5*(xn(nm,3)+xn(nm,1));yn(nm,6)=0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4)=0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4)=1; end
   efln(nm,5)=0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5)=1; end
   efln(nm,6)=0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6)=1; end

   nm=nm+1;
   xn(nm,1)=x(j,6); yn(nm,1)=y(j,6); efln(nm,1)=efl(j,6); %  third sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,3); yn(nm,3)=y(j,3); efln(nm,3)=efl(j,3);
   xn(nm,4)=0.5*(xn(nm,1)+xn(nm,2));yn(nm,4)=0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5)=0.5*(xn(nm,2)+xn(nm,3));yn(nm,5)=0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6)=0.5*(xn(nm,3)+xn(nm,1));yn(nm,6)=0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4)=0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4)=1; end
   efln(nm,5)=0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5)=1; end
   efln(nm,6)=0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6)=1; end

   nm=nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  fourth sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);
   xn(nm,4)=0.5*(xn(nm,1)+xn(nm,2));yn(nm,4)=0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5)=0.5*(xn(nm,2)+xn(nm,3));yn(nm,5)=0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6)=0.5*(xn(nm,3)+xn(nm,1));yn(nm,6)=0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4)=0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4)=1; end
   efln(nm,5)=0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5)=1; end
   efln(nm,6)=0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6)=1; end

 end % end of loop over current elements

 ne=4*ne;  % number of elements has increased by a factor of four
 for k=1:ne  % relabel the new points and put them in the master list                
    for l=1:6
     x(k,l)=xn(k,l); y(k,l)=yn(k,l);  efl(k,l)=efln(k,l);
    end
 end
end % end do of refinement loop
end % end if of refinement loop

% define the global nodes and the connectivity table
% six nodes of the first element are entered manually
p(1,1)=x(1,1); p(1,2)=y(1,1); gfl(1)=efl(1,1);
p(2,1)=x(1,2); p(2,2)=y(1,2); gfl(2)=efl(1,2);
p(3,1)=x(1,3); p(3,2)=y(1,3); gfl(3)=efl(1,3);
p(4,1)=x(1,4); p(4,2)=y(1,4); gfl(4)=efl(1,4);
p(5,1)=x(1,5); p(5,2)=y(1,5); gfl(5)=efl(1,5);
p(6,1)=x(1,6); p(6,2)=y(1,6); gfl(6)=efl(1,6);

c(1,1)=1;  % first  node of first element is global node 1
c(1,2)=2;  % second node of first element is global node 2
c(1,3)=3;  % third  node of first element is global node 3
c(1,4)=4;  % fourth node of first element is global node 4
c(1,5)=5;  % fifth  node of first element is global node 5
c(1,6)=6;  % sixth  node of first element is global node 6

ng=6;
% loop over further elements, Iflag=0 will signal a new global node
eps=0.000001;
for i=2:ne        % loop over elements
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
   ng=ng+1; p(ng,1)=x(i,j); p(ng,2)=y(i,j); gfl(ng)=efl(i,j);
   c(i,j)=ng;   % the jth local node of element is the new global node
 end
 end
end  % end of loop over elements
p(:,1)=lx*p(:,1)/nx; p(:,2)=ly*p(:,2); % scale the coordinates 
