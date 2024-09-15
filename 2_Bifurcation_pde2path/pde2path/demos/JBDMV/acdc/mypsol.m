function mypsol(dir,pt); 
p=loadp(dir,pt); y2=max(max(p.u(1:p.nu)),1); 
plotsol(p); axis([0 1 0 y2]); hold on; 
m=round(p.nu/2); p.u(m-2:m+4)', min(p.u(1:p.nu))
idx=find(p.u(1:p.nu)<1e-8); x=getpte(p); 
plot(x(idx),p.u(idx),'*m'); nola
