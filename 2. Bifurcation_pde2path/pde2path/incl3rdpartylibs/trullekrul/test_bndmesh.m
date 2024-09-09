clc; clear all; close all;
N = 10; [xy,tri,geomfunc] = mesh_rect(N); %options=gen_options(); bndmesh = gen_polygon_bndmesh(tri,xy,options);
[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri);
figure(); 	trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','o'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);

options = gen_options();
bndmesh = bndmesh_polygon(tri,xy,[],options);
% PLOT BNDMESH
hold on; for i=1:size(bndmesh.edg,1)
    text((xy(bndmesh.edg(i,1),1)+xy(bndmesh.edg(i,2),1))/2, (xy(bndmesh.edg(i,1),2)+xy(bndmesh.edg(i,2),2))/2,num2str(bndmesh.IDs(i)),'color','r','FontWeight','bold','fontsize',14);end; hold off;

%% PLOT BNDMESH
%nbndedg = bks_bndedg2edg(edg,nd2edg,bndmesh);
%hold on; for i=1:size(nbndedg,1)
    %text((xy(edg(nbndedg(i),1),1)+xy(edg(nbndedg(i),2),1))/2, (xy(edg(nbndedg(i),1),2)+xy(edg(nbndedg(i),2),2))/2,num2str(bndmesh.IDs(i)),'color','r','FontWeight','bold','fontsize',14);end; hold off;
    %
%figure();
%hold on; for i=1:size(bndmesh.fac,1)
    %text((xy(bndmesh.fac(i,1),1)+xy(bndmesh.fac(i,2),1)+xy(bndmesh.fac(i,2),3))/3, (xy(bndmesh.fac(i,1),2)+xy(bndmesh.fac(i,2),2)+xy(bndmesh.fac(i,2),3))/3,num2str(bndmesh.IDs(i)),'color','r','FontWeight','bold','fontsize',14);end; hold off;