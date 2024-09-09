clc; clear all; close all;
N = 15; test_box_uniform2D;
options.qualM = 1;
triQ = elem_qual(tri,xy,Nmetric,options);
tri_ = reshape(1:3*size(tri,1),3,size(tri,1))'; xy_ = xy(reshape(tri',numel(tri),1),:);
figure(2); trisurf(tri_,xy_(:,1),xy_(:,2),reshape(repmat(triQ',3,1),3*size(tri,1),1),'edgecolor','none');
set(gca,'CameraTarget',[0.5,0.5,0.]); set(gca,'Camerapositio',[0.5,0.5,1.]); colorbar(); axis off; box off;
set(gca,'fontsize',16); title('2D Vassilevski Functional');
%print('-dpng','elem_qual.png'); print('-depsc2','elem_qual.eps');
