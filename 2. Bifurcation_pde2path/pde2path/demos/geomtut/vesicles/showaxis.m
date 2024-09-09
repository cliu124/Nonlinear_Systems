function showaxis(p) % show the axis used as rotational axis, see rotax.m 
A=4*[p.w;p.om;p.rh]; pplot(p,10); figure(10);  
quiver3(zeros(3,1),zeros(3,1),zeros(3,1), A(:,1), A(:,2),A(:,3),'Color','r'); 