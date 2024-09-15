clc; clear all; close all;
N = 10; [xy,tri,geomfunc] = mesh_rect(N); 
[edg,edg2tri,tri2edg,nd2tri,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri); 
	trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);

ii = [1 2]; % [1 2 3];

if any(ii==1)
 % NODE NUMBERS
 hold on; for i=1:size(xy,1)
     text(xy(i,1),xy(i,2),num2str(i),'color','b','FontWeight','bold','fontsize',14); end; hold off;
end;
if any(ii==2)
 % EDGE NUMBERS
 hold on; for i=1:size(edg,1)
     text((xy(edg(i,1),1)+xy(edg(i,2),1))/2,(xy(edg(i,1),2)+xy(edg(i,2),2))/2,num2str(i),'color','r','FontWeight','bold','fontsize',14); end; hold off;
end;
if any(ii==3)
 % ELEMENT NUMBERS
 hold on; for i=1:size(tri,1)
     text((xy(tri(i,1),1)+xy(tri(i,2),1)+xy(tri(i,3),1))/3,(xy(tri(i,1),2)+xy(tri(i,2),2)+xy(tri(i,3),2))/3,num2str(i),'color','k','FontWeight','bold','fontsize',14); end; hold off;
end;