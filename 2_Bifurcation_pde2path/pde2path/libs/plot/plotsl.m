function plotsl(o,u,step) 
% PLOTSL: plot slice of u through grid.
%
%  plotsl(o,u,step)
zslice= unique(o.p(3,:)); 
bdtri = find(o.e(4,:)==0); 
bdtri = bdtri(1:length(bdtri)/2);
                       
mx = max(u); mm = min(u);    
for k = 1:step:length(zslice)
   indx=o.p(3,:)==zslice(k); % find points in slice
   cf = u(indx); % fuval in slice 
   X1 = zeros(3,length(bdtri));
   X2 = zeros(3,length(bdtri)); 
   cind = zeros(1,length(bdtri));
   for k2 = bdtri
     c = sum(cf(o.e(1:3,k2)))/3; % average over edges 
     cind(k2)=c; 
     X1(:,k2)=o.p(1,o.e(1:3,k2))'; X2(:,k2)=o.p(2,o.e(1:3,k2))';
   end 
   h = patch(X1,zslice(k)*ones(3,length(bdtri)),X2,cind);set(h,'EdgeColor','none');               
end
axis equal; view(55,30)
end
