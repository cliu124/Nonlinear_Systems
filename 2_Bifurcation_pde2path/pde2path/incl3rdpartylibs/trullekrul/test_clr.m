clc; clear all; close all;
N = 15; test_box_uniform2D; %N = 10; [xy,tri,geomfunc] = mesh_rect(N); 
[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri);
options = gen_options();
options.greedyCLR = -1;
ii = 1;
trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);
axis off; box off;
thecolors = [0.0000    0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840
             0.0000    0.0000    0.0000];
%ocolors = jet; thecolors = [thecolors; ocolors(10:12:60,:)];

switch ii
case 0
	edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
	%edg2clr = gen_edgclr(edg_ngh,options);
	edg2clr = bks_clr(edg_ngh,options);
	hold on; for i=1:size(edg,1)
    text((xy(edg(i,1),1)+xy(edg(i,2),1))/2,(xy(edg(i,1),2)+xy(edg(i,2),2))/2,num2str(edg2clr(i)),'color','b','FontWeight','bold','fontsize',14); end; hold off;
case 1
	%ngh = bks_nd2tri2nd(edg,nd2edg);
	ngh = bks_nd2tri2nd2(tri,nd2tri);
	nd2clr = bks_clr(ngh,options);
	 %hold on; for i=1:size(xy,1)
     %text(xy(i,1),xy(i,2),num2str(nd2clr(i)),'color','r','FontWeight','bold','fontsize',14); end; hold off;
     hold on; for i=1:max(nd2clr)
     plot(xy(i==nd2clr,1),xy(i==nd2clr,2),'o','markerfacecolor',thecolors(i,:),'markeredgecolor',thecolors(i,:)); end; hold off;
     %print('-dpng','box_color.png'); print('-depsc2','box_color.eps');
case 2
	[nghOa,ngh] = bks_edg2nd2tri2nd2edg(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,options);
		%edg2clr = gen_edgclr(ngh,options);
	edg2clr = bks_clr(ngh,options);
	hold on; for i=1:size(edg,1)
     text((xy(edg(i,1),1)+xy(edg(i,2),1))/2,(xy(edg(i,1),2)+xy(edg(i,2),2))/2,num2str(edg2clr(i)),'color','b','FontWeight','bold','fontsize',14); end; hold off;
end;
