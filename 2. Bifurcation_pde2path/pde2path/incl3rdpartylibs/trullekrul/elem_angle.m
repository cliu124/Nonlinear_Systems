function [angle,badedg,badangle] = elem_angle(tri,xy,options)
v1 = xy(tri(:,2),:) - xy(tri(:,1),:);
v2 = xy(tri(:,3),:) - xy(tri(:,1),:);
v3 = xy(tri(:,2),:) - xy(tri(:,3),:);
if size(tri,2) == 3
	a = [-0.03897865137325461082351552022373653016984462738037109375  0.4217895987591873119271212999592535197734832763671875  -1.667050099615339231462485258816741406917572021484375  2.28287871497437411250075456337071955204010009765625];
	%a = [-0.03897865 0.4217896 -1.66705 2.2828787];

	v1 = v1./repmat(sqrt(sum(v1.^2,2)),1,size(xy,2));
	v2 = v2./repmat(sqrt(sum(v2.^2,2)),1,size(xy,2));
	v3 = v3./repmat(sqrt(sum(v3.^2,2)),1,size(xy,2));
	angle1 = acos(sum(v1.*v2,2));
	angle2 = acos(sum(v1.*v3,2));
	%if sum(abs(imag(angle1))) ~= 0 || sum(abs(imag(angle2))) ~= 0
	%	error('in angle calculation');
	%end;
	angle3 = pi-angle1-angle2;
	%mean([angle1 angle2 angle3])
	%myf = @(x_) a(1)*x_.^4 + a(2)*x_.^3 + a(3)*x_.^2 + a(4)*x_;	
	%angle = myf(max([angle1 angle2 angle3],[],2));
	if options.qualP == eps && nargout == 1
		angle = [angle1; angle2; angle3];
		return;
	end;
	if nargout == 3
		I = [angle1 angle2 angle3] > options.qualP;
		badedg = sort([tri(I(:,1),3) tri(I(:,1),2); ...
 	  		       tri(I(:,2),1) tri(I(:,2),3); ...
 	  		       tri(I(:,3),2) tri(I(:,3),1)],2);
 	        badangle = [angle1(I(:,1)); angle2(I(:,2)); angle3(I(:,3))];
 	end;

	myf = @(x_) min([ones(size(x_)) -1/(pi-options.qualP)*(x_-options.qualP)+1],[],2);
	angle = myf([angle1 angle2 angle3]);
else
	v4 = xy(tri(:,1),:)-xy(tri(:,4),:);
	v5 = xy(tri(:,2),:)-xy(tri(:,4),:);
	%v6 = xy(tri(:,3),:)-xy(tri(:,4),:);
	a = [-0.01249403951185311910376807276179533801041543483734130859375 0.1928840463398532600880486143068992532789707183837890625 -1.0780219448064494169869931283756159245967864990234375 1.8704102355158995774075947338133119046688079833984375];
	%a = [-0.012494 0.192884 -1.078022 1.870410];
	v123 = cross(v1,v2,2); v123 = v123./repmat(sqrt(sum(v123.^2,2)),1,3);
	v124 = cross(v1,v4,2); v124 = v124./repmat(sqrt(sum(v124.^2,2)),1,3);
	v134 = cross(v2,v4,2); v134 = v134./repmat(sqrt(sum(v134.^2,2)),1,3); 
	v234 = cross(v3,v5,2); v234 = v234./repmat(sqrt(sum(v234.^2,2)),1,3);
	angle1 =  acos(sum(-v123.*v124,2));
	angle2 =  acos(sum(v123.*v134,2));
	angle3 =  acos(sum(v123.*v234,2));
	angle4 =  acos(sum(v124.*v134,2));
	angle5 =  acos(sum(v124.*v234,2));
	angle6 =  acos(sum(-v134.*v234,2));
	%[angle1 angle2 angle3 angle4 angle5 angle6]*180/pi
	%mean([angle1 angle2 angle3 angle4 angle5 angle6]*180/pi)
	%myf = @(x_) a(1)*x_.^4 + a(2)*x_.^3 + a(3)*x_.^2 + a(4)*x_;
	%angle = myf(max([angle1 angle2 angle3 angle4 angle5 angle6],[],2));
	if options.qualP == eps && nargout == 1
		angle = [angle1; angle2; angle3; angle4; angle5; angle6];
		return;
	end;
	if nargout == 3
		I = [angle1 angle2 angle3 angle4 angle5 angle6] > options.qualP;
		badedg = sort([tri(I(:,1),3) tri(I(:,1),4); ...
 	  		       tri(I(:,2),4) tri(I(:,2),2); ...
 	  		       tri(I(:,3),1) tri(I(:,3),4); ...
 	  		       tri(I(:,4),2) tri(I(:,4),3); ...
 	  		       tri(I(:,5),3) tri(I(:,5),1); ...
 	  		       tri(I(:,6),1) tri(I(:,6),2)],2);
 	        badangle = [angle1(I(:,1)); angle2(I(:,2)); angle3(I(:,3)); ...
 	                    angle4(I(:,4)); angle5(I(:,5)); angle6(I(:,6))];
	end;
	myf = @(x_) min([ones(size(x_)) -1/(pi-options.qualP)*(x_-options.qualP)+1],[],2);
	angle = myf([angle1 angle2 angle3 angle4 angle5 angle6]);
end;



%a(1)*x^3+a(2)*x^2+a(3)*x+a(4)
%a(4) = 0
%3*a(1)*x1^2+2*a(2)*x1+a(3) = 0
%a(1)*x1^3+a(2)*x1^2+a(3)*x1 = 1
%a(1)*pi^3+a(2)*pi^2+a(3)*pi = 0
%x1 = pi/3; A=[3*x1^2 2*x1 1; x1^3 x1^2 x1; pi^3 pi^2 pi]; b = [0;1;0]; sprintf('%1.64f %1.64f %1.64f',A\b)
%x1 = acos(1/3); A=[3*x1^2 2*x1 1; x1^3 x1^2 x1; pi^3 pi^2 pi]; b = [0;1;0]; sprintf('%1.64f %1.64f %1.64f',A\b)
%a = A\b; myf = @(x_) a(1)*x_.^3 + a(2)*x_.^2 + a(3)*x_; x = linspace(0,pi,101); plot(x,myf(x))

%a(0)*x^4+a(1)*x^3+a(2)*x^2+a(3)*x+a(4)
%a(4) = 0
%4*x(0)*x1^3+3*a(1)*x1^2+2*a(2)*x1+a(3) = 0
%a(0)*x1^4+a(1)*x1^3+a(2)*x1^2+a(3)*x1 = 1
%a(0)*pi^4+a(1)*pi^3+a(2)*pi^2+a(3)*pi = 0
%12*a(0)*pi^2+6*a(1)*pi+2*a(2)= 0
%x1 = pi/3; A=[4*x1^3 3*x1^2 2*x1 1; x1^4 x1^3 x1^2 x1; pi^4 pi^3 pi^2 pi; 12*pi^2 6*pi 2 0]; b = [0;1;0; 0]; sprintf('%1.64f %1.64f %1.64f %1.64f',A\b)
%x1 = acos(1/3); A=[4*x1^3 3*x1^2 2*x1 1; x1^4 x1^3 x1^2 x1; pi^4 pi^3 pi^2 pi; 12*pi^2 6*pi 2 0]; b = [0;1;0; 0]; sprintf('%1.64f %1.64f %1.64f %1.64f',A\b)
%a = A\b; myf = @(x_) a(1)*x_.^4 + a(2)*x_.^3 + a(3)*x_.^2 + a(4)*x_; x = linspace(0,pi,101); plot(x,myf(x))