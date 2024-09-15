function i=triint(u,p,t)
% TRIINT: integrate u over Om given by p(oints) and t(riangles), obsolete 
%
%  i=triint(u,p,t)
%  u=nodal values (scalar, i.e., size(u)=np x 1) 
ndim=size(p,1); 
if ndim==1; 
   n=size(p,2); um=0.5*(u(1:n-1)+u(2:n)); dx=p(2:n)-p(1:n-1); 
   i=sum(dx'.*um); 
else if ndim==2;
   a1=t(1,:);a2=t(2,:);a3=t(3,:); % Corner point indices 
   % Triangle sides
   r23x=p(1,a3)-p(1,a2); r23y=p(2,a3)-p(2,a2);
   r31x=p(1,a1)-p(1,a3); r31y=p(2,a1)-p(2,a3);
   ar=abs(r31x.*r23y-r31y.*r23x)/2;% areas of tiangles 
   ut=pdeintrp(p,t,u); i=ut*ar';
   end
end

