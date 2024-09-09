function [c,h]=levplot(dir,pt,cl) 
% levplot: plot level-lines u1=cl  
p=loadp(dir,pt); gr=p.pdeo.grid; x=gr.p(1,:); y=gr.p(2,:); t=gr.t(1:3,:); t=t'; 
[c,h]=tricontour(t,x,y,p.u(1:p.np),[cl,cl]); hold on 
for k=1:length(gr.e(5,:))
    line(gr.p(1,gr.e(1:2,k)),gr.p(2,gr.e(1:2,k)), 'LineWidth',1,'color','blue')
end 
if size(cl,2)==2
if cl(1)==cl(2); tit=[dir '/' pt ',lev=' mat2str(cl(1),3)];
else tit=[dir '/' pt ',lev=(' mat2str(cl(1),2) ',' mat2str(cl(2),2) ')']; 
end 
title(tit); 
else clabel(c); 
end 