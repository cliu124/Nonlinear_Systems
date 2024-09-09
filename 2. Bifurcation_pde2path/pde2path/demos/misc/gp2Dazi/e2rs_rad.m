function [p,idx]=e2rs_rad(p,u) % elements to refine selector for refinement near (x,y)=0 (r=0) 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; 
p2c=point2CenterMatrix(p.pdeo.grid);
xt=p2c*x; yt=p2c*y; 
r=xt.^2+yt.^2; 
[rs,idx1]=sort(r); idx=idx1(1:p.nref); 