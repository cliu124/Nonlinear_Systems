function r=myfilter(er)
% myfilter: helper function for ampsys 
% put equal rows of er(1:n,1:m-1) together and sum entries of er(:,m) 
n=size(er,1); m=size(er,2); r=[];
for i=1:n
   id=[];
   for j=1:n
       if er(i,1:m-1)==er(j,1:m-1); id=[id,j]; end
   end
   r=[r;[er(i,1:m-1),sum(er(id,m))]];
end
r=unique(r,'rows'); 

