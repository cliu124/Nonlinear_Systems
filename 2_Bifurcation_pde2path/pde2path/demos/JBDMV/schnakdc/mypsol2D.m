function mypsol2D(dir,pt,v); 
p=loadp(dir,pt); y2=max(max(p.u(1:p.nu)),1); 
plotsol(p); colormap parula; hold on; 
m=round(p.nu/2); p.u(m-2:m+4)', min(p.u(1:p.nu))
idx=find(p.u(1:p.nu)<1e-8); x=getpte(p); 
plot3(x(1,idx),x(2,idx),p.u(idx),'*m'); 
nola; view(v); getlam(p)