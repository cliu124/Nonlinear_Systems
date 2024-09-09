function mmI = metric_avg(edg,Nmetric,options);
if options.log_m_add
	if size(edg,2) == 2
		mmI = add_log_metric(Nmetric(edg(:,1),:),Nmetric(edg(:,2),:));
	elseif size(edg,2) == 3
		mmI = 	add_log_metric(Nmetric(edg(:,1),:),Nmetric(edg(:,2),:),Nmetric(edg(:,3),:));
	else
		mmI = 	add_log_metric(Nmetric(edg(:,1),:),Nmetric(edg(:,2),:),Nmetric(edg(:,3),:),Nmetric(edg(:,4),:));
	end;
else
	if size(edg,2) == 2
		mmI = (Nmetric(edg(:,1),:)+Nmetric(edg(:,2),:))/2;
	elseif size(edg,2) == 3
		mmI = (Nmetric(edg(:,1),:)+Nmetric(edg(:,2),:)+Nmetric(edg(:,3),:))/3;
	else
		mmI = (Nmetric(edg(:,1),:)+Nmetric(edg(:,2),:)+Nmetric(edg(:,3),:)+Nmetric(edg(:,4),:))/4;
	end;
end;


function out_sum = add_log_metric(metric1,metric2,metric3,metric4)
%[lambda1 ,lambda2 ,v1xn ,v1yn ] = analyt_eig(metric1);
%mm_log1 = [v1xn.^2.*log(lambda1)+v1yn.^2.*log(lambda2) ...
       %v1yn.*v1xn.*(log(lambda1)-log(lambda2)) ...
       %v1yn.^2.*log(lambda1)+v1xn.^2.*log(lambda2)];
%[lambda1,lambda2,v1xn,v1yn] = analyt_eig(metric2);
%mm_log2 = [v1xn.^2.*log(lambda1)+v1yn.^2.*log(lambda2) ...
       %v1yn.*v1xn.*(log(lambda1)-log(lambda2)) ...
       %v1yn.^2.*log(lambda1)+v1xn.^2.*log(lambda2)];
%logsum = mm_log1+mm_log2;
%if nargin == 3
	%[lambda1,lambda2,v1xn,v1yn] = analyt_eig(metric3);
	%logsum = logsum+[v1xn.^2.*log(lambda1)+v1yn.^2.*log(lambda2) ...
	       %v1yn.*v1xn.*(log(lambda1)-log(lambda2)) ...
	       %v1yn.^2.*log(lambda1)+v1xn.^2.*log(lambda2)];
%end;
%logsum = logsum/nargin;
%[lambda1,lambda2,v1xn,v1yn] = analyt_eig(logsum);
%out_sum = [v1xn.^2.*exp(lambda1)+v1yn.^2.*exp(lambda2) ...
       %v1yn.*v1xn.*(exp(lambda1)-exp(lambda2)) ...
       %v1yn.^2.*exp(lambda1)+v1xn.^2.*exp(lambda2)];
       
[eigL, eigR] = analyt_eig(metric1);
mm_log1 = analyt_prod(analyt_fulleig(log(eigL)),eigR);
[eigL, eigR] = analyt_eig(metric2);
mm_log2 = analyt_prod(analyt_fulleig(log(eigL)),eigR);
logsum = mm_log1+mm_log2;
if nargin == 3
	[eigL, eigR] = analyt_eig(metric3);
	logsum = logsum+analyt_prod(analyt_fulleig(log(eigL)),eigR);
end;
if nargin == 4
	[eigL, eigR] = analyt_eig(metric4);
	logsum = logsum+analyt_prod(analyt_fulleig(log(eigL)),eigR);
end;
logsum = logsum/nargin;
[eigL, eigR] = analyt_eig(logsum);
out_sum = analyt_prod(analyt_fulleig(exp(eigL)),eigR);