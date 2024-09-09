function [badtri,z] = elem_inv(tri,xy,xy1,xy2,xy3,xy4)
if (nargin == 2 && numel(tri) == 0) || (nargin ~= 2 && numel(xy1) == 0) 
    badtri = []; z = []; 
    return;
end;
if nargin == 3 && size(tri,2) == 3 && size(xy,2) == 3
	nvec = xy1;
end;
% orient for check of inversion later on:
%xy=xy'; 
if nargin <= 3
	xy1 = xy(tri(:,1),:);
	xy2 = xy(tri(:,2),:);
	xy3 = xy(tri(:,3),:);
	if size(tri,2) == 4
		xy4 = xy(tri(:,4),:);
	end;
end;
v1 = xy2-xy1;
v2 = xy3-xy1;
if size(xy1,2) == 2
	z = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
else
	crosspro = cross(v1,v2,2); %[v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2) ...
%			   v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3) ...
%			   v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)];
	if (size(tri,2) == 3 || nargin == 5) && size(xy1,2) == 3 
		if nargin == 3
		z = sum(nvec.*crosspro,2); %exterior face taking normal into account
		else
		z = sqrt(sum(crosspro.^2,2)); %interior face only requires area calculation
		end;
	else
		v3 = xy4 - (xy1+xy2+xy3)/3;
	         z = sum(v3.*crosspro,2);
        end;
end;
badtri = z<0; %counter-clock-wise