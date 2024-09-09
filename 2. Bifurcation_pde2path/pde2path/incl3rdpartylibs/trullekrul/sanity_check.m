function [] = sanity_check(tri,xy,triQ,Nmetric,options)
[I,area] = elem_inv(tri,xy);
if any(I)
    trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]); %save('for_debug.mat');
    error(sprintf('inverted element (%1.1f,%1.1f,%1.1f)',mean(xy(tri(find(I,1),:),:))));
end;
if nargin == 2
	return
end;
if options.debug == 2
%if options.qualM == 9
%triQ_ = elem_qual(tri, xy, metric_iso(xy, 1), options);
%else
triQ_ = elem_qual(tri, xy, Nmetric,options);
%end;
if any(abs(triQ_-triQ)>1e-8)
	error(sprintf('flawed quality list (Emax=%0.0e)',max(abs(triQ_-triQ))));
end;
end;
if size(xy,2) == 3
	area = area/6;
else
	area = area/2;
end;
if options.debug == 2 && options.area ~= 0 && abs(options.area-sum(area))>max(options.minA,1e-11)
    error_ = abs(options.area-sum(abs(area)));
    if size(xy,2) == 2
    trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]); %save('for_debug.mat');
    x = xy(:,1); y = xy(:,2);
%     hold on; for i=1:size(tri,1)
%         text((x(tri(i,1))+x(tri(i,2))+x(tri(i,3)))/3,(y(tri(i,1))+y(tri(i,2))+y(tri(i,3)))/3,num2str(i),'color','r','FontWeight','bold','fontsize',14); end; hold off;    
    hold on; plot(mean(x(tri),2),mean(y(tri),2),'.r'); hold off; end;
    error(sprintf('domain area reduced (%0.0e), (curved geometry? set options.area=0)',error_));
end;