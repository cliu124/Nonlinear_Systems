function p=tripledrop(p)
% tripledrop: triple the size of p.mat.drop and p.mat.fill  
r=size(p.mat.drop,1); s=size(p.mat.drop,2);
if r>1;  zrs=sparse(r,s); drop=p.mat.drop; fill=p.mat.fill; 
  p.mat.drop=[[drop zrs zrs]; [zrs drop zrs]; [zrs zrs drop]]; 
  p.mat.fill = [[fill zrs' zrs']; [zrs' fill zrs']; [zrs' zrs' fill]]; 
end