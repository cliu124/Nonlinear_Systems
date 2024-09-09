j=2; figure(j); clf;  s2=sqrt(2); s3=sqrt(3); %s2=1.6;
x=-pi:pi/10:pi; if j==1; x=s3*x; else; x=s2*x; end  
y=x; z=x; [X,Y,Z]=meshgrid(x,y,z);  
lev=[-0.5 0.5]; col=['b' 'r']; s2=sqrt(2); s3=sqrt(3); %s2=1.6;
switch j
    case 1; v=cos((X+Y+Z)/s3)+cos((X-Y-Z)/s3)+cos((-X+Y-Z)/s3)+cos((-X-Y+Z)/s3); % FCC
    case 2; v=1*cos((X+Y)/s2)+0*cos((Y+Z)/s2)+0*cos((X+Z)/s2)+...
            1*cos((X-Y)/s2)+0*cos((Y-Z)/s2)+0*cos((-X+Z)/s2); 
end 
%v=cos(X)+cos(Y)+cos(Z); % SC
for i=1:2; 
  ip=patch(isosurface(X,Y,Z,v,lev(i))); isonormals(X,Y,Z,v,ip); 
 set(ip,'FaceColor',col(i),'EdgeColor','none','facealpha',0.8); hold on;
end 
%daspect([1,1,1]); 
view([-40 60]); axis tight
camlight 
lighting phong