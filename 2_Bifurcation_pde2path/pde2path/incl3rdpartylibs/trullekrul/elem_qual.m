function [quality,c] = elem_qual(tri,xy,Nmetric,options,nvec)
%options.qualM == 1 Vassilevski
%options.qualM == 2 Orth, LW(H), worst
%options.qualM == 3 Orth, LW(H), worst vassilevski
%options.qualM == 4 Orth, ellipse, worst
%options.qualM == 5 Orth, ellipse, worst vassilevski
%options.qualM == 6 Orth, LW(H), alt worst
%options.qualM == 7 Orth, ellipse, alt worst
%options.qualM == 8 Home brew condition (sliver functional), now it is shephaoprd
%options.qualM == 9 Condition number

if nargin == 5 && size(tri,2) ~= 2 %we wont bother calculating quality of inverted elements
  quality = repmat(-1,size(tri,1),1);
  if size(nvec,2) == 3
    [Ibad,area] = elem_inv(tri,xy,nvec); 
  else
    [Ibad,area] = elem_inv(tri,xy); 
  end;
  if size(tri,2) == 4 || size(tri,2) == 3
  	Ibad = area < options.minA;
  end;
  if any(not(Ibad))
  quality(not(Ibad)) = elem_qual_slow(tri(not(Ibad),:),xy,Nmetric,options);
  end;
else
  if nargout == 2
	[quality,c] = elem_qual_slow(tri,xy,Nmetric,options);
  else
  	quality = elem_qual_slow(tri,xy,Nmetric,options);
  end;
end;
	
function [quality,c] = elem_qual_slow(tri,xy,Nmetric,options)	
if size(tri,2) == 2
	if nargout == 2
[quality,c] = gen_rel_edgL(tri,xy,Nmetric,options);
	else
	quality = gen_rel_edgL(tri,xy,Nmetric,options);
	end;
else
        if options.qualP > 0
  	  quality = elem_angle(tri,xy,options);
  	  return;
  	end;
if options.qualM == 1 || options.qualM > 7 
quality = vassilevski(tri,xy,Nmetric,options);
end;
if options.qualM == 2 || options.qualM == 3 || options.qualM == 6
    quality = worst_HL(tri,xy,Nmetric,options);
end;
if options.qualM == 4 || options.qualM == 5 || options.qualM == 7
     quality = steiner_ell_metric(tri,xy,Nmetric,options);
end;
end;


function quality = vassilevski(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum((v1-v2).^2,2)]);
else
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum(v3.^2,2) sum((v1-v2).^2,2) sum((v1-v3).^2,2) sum((v2-v3).^2,2)]);
end;
if options.qualM == 10
quality = area; return;
end;
if size(tri,2) == 3
if options.qualM == 1
quality = area./mean(L,2).^2.*myf(mean(L,2)); 
elseif options.qualM == 8
%quality = area./mean(L,2).^2;
quality = area./mean(L.^2,2);
else %options.qualM == 9
quality = area./mean(L.^2,2);
end;
else %3D
	[areaA1,H1] = gen_crosspro(xy1I,xy2I,xy3I,xy4I);
	[areaA2,H2] = gen_crosspro(xy4I,xy2I,xy3I,xy1I);
	[areaA3,H3] = gen_crosspro(xy1I,xy4I,xy3I,xy2I);
        [areaA4,H4] = gen_crosspro(xy1I,xy2I,xy4I,xy3I);
        L2 = [areaA1 areaA2 areaA3 areaA3];
        if options.qualM == 1
        %quality = area./mean(L,2).^3.*myf(mean(L2,2));
        quality = area./mean(L,2).^3.*myf(mean(L,2));
        %quality = area./mean(L.^2,2).^(3/2);
        %quality = area./mean(L2,2).^(3/2).*myf(mean(L,2));
        %quality = area./mean(L2,2).^(3/2).*mean(myf(L),2);
        elseif options.qualM == 8
       	%quality = area./mean(L,2).^3; 
       	%quality = area./mean(L,2)./mean(L2,2); 
       	quality = area.^2./mean(L.^2,2).^3; 
        else %options.qualM == 9
        %quality = area./mean(L,2).^3.*myf(mean(L2,2)); 
        quality = area./sqrt(mean(L.^2,2).*mean(L2.^2,2));
        %quality = area./mean(L.^2,2).^(3/2); %frenchy
        end;
	%qual2Ds = elem_qual([tri(:,1:3); tri(:,2:4); tri(:,[1 3 4]); tri(:,[1 2 4])],xy,Nmetric,options); qual2Ds = reshape(qual2Ds,size(tri,1),4);
	%quality = area./mean(L2,2).^(3/2).*mean(qual2Ds,2);
end;


function quality = steiner_ell_metric(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2] = tri_in_metric_space(tri,xy,Nmetric,options);
	XI = elem_metric(tri,xy,xy1I,xy2I,xy3I);

else
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options);
	XI = elem_metric([],[],xy1I,xy2I,xy3I,xy4I);
end;
[eigL,eigR] = analyt_eig(XI); 
if any(eigL(:)<0)
    warning('neg eig'); eigL = abs(eigL);
end;
eigL = 1./sqrt(eigL);
if options.qualM == 4
	quality = min(myfa(eigL),[],2).*sign(area);
elseif options.qualM == 7
	quality = min(relE(eigL),[],2).*sign(area);
else
	quality = area./mean(eigL,2).^(size(tri,2)-1).*myf(mean(eigL,2));
end;


function quality = worst_HL(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum((v1-v2).^2,2)]);	
	H = repmat(abs(area),1,3)./L;
	%2D unit element has height sqrt(3)/2
	vals = [L(:,1) H(:,1) L(:,2) H(:,2) L(:,3) H(:,3)];
else
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum(v3.^2,2) sum((v1-v2).^2,2) sum((v1-v3).^2,2) sum((v2-v3).^2,2)]); %(12,13,14,23,24,34)
	[areaA1,H1] = gen_crosspro(xy1I,xy2I,xy3I,xy4I);
	[areaA2,H2] = gen_crosspro(xy4I,xy2I,xy3I,xy1I);
	[areaA3,H3] = gen_crosspro(xy1I,xy4I,xy3I,xy2I);
        [areaA4,H4] = gen_crosspro(xy1I,xy2I,xy4I,xy3I);
	Wa = [areaA1 areaA1 areaA3 areaA1 areaA2 areaA2]./L;
	Wb = [areaA4 areaA3 areaA4 areaA2 areaA4 areaA3]./L;
	Ha = [H1 H1 H3 H1 H2 H2]; %Ha = repmat(area,1,6)./(L.*Wa);
	Hb = [H4 H3 H4 H2 H4 H3]; %Hb = repmat(area,1,6)./(L.*Wb);
	vals = [L Wa Wb Ha Hb];
end;
if options.qualM == 2 || options.qualM == 6
	nzero = area~=0;
	if options.qualM == 2 
	worst_ = min(myfa(vals(nzero,:)),[],2);
	else
	worst_ = min(relE(vals(nzero,:)),[],2);
	end;
	worst = zeros(size(L,1),1); worst(nzero) = worst_;
	quality = worst.*sign(area);
else
	myf2 = @(x_) area./mean(x_,2).^(size(tri,2)-1).*myf(mean(x_,2));
	if size(tri,2) == 3
		qual = [myf2([L(:,1) H(:,1)]) myf2([L(:,2) H(:,2)]) myf2([L(:,3) H(:,3)])];
	else
		qual = [myf2([L(:,1) Wa(:,1) Ha(:,1)]) myf2([L(:,1) Wb(:,1) Hb(:,1)]) ...
         		        myf2([L(:,2) Wa(:,2) Ha(:,2)]) myf2([L(:,2) Wb(:,2) Hb(:,2)]) ...
         		        myf2([L(:,3) Wa(:,3) Ha(:,3)]) myf2([L(:,3) Wb(:,3) Hb(:,3)]) ...
         		        myf2([L(:,4) Wa(:,4) Ha(:,4)]) myf2([L(:,4) Wb(:,4) Hb(:,4)]) ...         		        
         		        myf2([L(:,5) Wa(:,5) Ha(:,5)]) myf2([L(:,5) Wb(:,5) Hb(:,5)]) ...
         		        myf2([L(:,6) Wa(:,6) Ha(:,6)]) myf2([L(:,6) Wb(:,6) Hb(:,6)])];
	end;
	quality = min(qual,[],2);%.*sign(area);
end;

function [mm,xy1,xy2,xy3,xy4] = calc_tri_metric(tri,xy,Nmetric,options);
xy1 = xy(tri(:,1),:); xy2 = xy(tri(:,2),:); xy3 = xy(tri(:,3),:); 
if size(tri,2) == 4
	xy4 = xy(tri(:,4),:);
end;
mm = metric_avg(tri,Nmetric,options);

function out = myf(x)
out = (min(x,1./x).*(2-min(x,1./x))).^3;

function out = myfa(myx)
out= min(myx,1./myx);

function out = relE(myx)
out = (myx.*(myx<1)+(2-myx).*(myx>=1));

function [crosspro,H] = gen_crosspro(x1,x2,x3,x4)
v1 = x2-x1;
v2 = x3-x1;
v3 = x4 - (x1+x2+x3)/3.;
H = sqrt(sum(v3.^2,2))/sqrt(2/3); %3D unit element has height sqrt(2/3)
crosspro = cross(v1,v2,2); %[v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2) ...
%		   v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3) ...
%		   v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)];
crosspro = sqrt(sum(crosspro.^2,2))*0.5/(sqrt(3)/4);

function [edgnowl,c] = gen_rel_edgL(edg,xy,Nmetric,options)
mmI = metric_avg(edg,Nmetric,options);
x = xy(:,1); y = xy(:,2);
x_ = reshape(x(edg),size(edg)); y_ = reshape(y(edg),size(edg));
	
if size(xy,2) == 3 %3D
	z = xy(:,3); z_ = reshape(z(edg),size(edg));
	xy1_ = [x_(:,1) y_(:,1) z_(:,1)]; xy2_ = [x_(:,2) y_(:,2) z_(:,2)];
else
	xy1_ = [x_(:,1) y_(:,1)]; xy2_ = [x_(:,2) y_(:,2)];
end;

%xyI1 = analyt_prod(xy1_,mmI); 
%xyI2 = analyt_prod(xy2_,mmI);
%t = xyI1-xyI2;
%edgnowl = sqrt(sum(t.^2,2)); %edge length
t = xy1_-xy2_;
edgnowl = sqrt(sum(analyt_prod(t,mmI).^2,2)); %edge length

if nargout == 2
	if options.spltc == 1
		c = (xy1_ + xy2_)/2;
	else
	%t = xy1_-xy2_; %t = t./repmat(sqrt(sum(t.^2,2)),1,size(t,2)); %compute invLs
	mmI1 = Nmetric(edg(:,1),:); L1 = sqrt(sum(analyt_prod(t,mmI1).^2,2));
	mmI2 = Nmetric(edg(:,2),:); L2 = sqrt(sum(analyt_prod(t,mmI2).^2,2));
	c = (xy1_.*repmat(L1,1,size(xy1_,2))+xy2_.*repmat(L2,1,size(xy1_,2)))./repmat(L1+L2,1,size(xy1_,2));
	end;
end;

function [Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[mm,xy1,xy2,xy3] = calc_tri_metric(tri,xy,Nmetric,options);
else
	[mm,xy1,xy2,xy3,xy4] = calc_tri_metric(tri,xy,Nmetric,options);
end;
xy1I = analyt_prod(xy1,mm);
xy2I = analyt_prod(xy2,mm);
xy3I = analyt_prod(xy3,mm);
v1 = xy2I - xy1I; 
v2 = xy3I - xy1I;
if size(tri,2) == 4
	xy4I = analyt_prod(xy4,mm);
	v3 = xy4I - xy1I;
	[Ibad,area] = elem_inv([],[],xy1I,xy2I,xy3I,xy4I);
	area = (1/6)/(sqrt(2)/12)*area; %volume
else
	[Ibad,area] = elem_inv([],[],xy1I,xy2I,xy3I);
	area = 0.5/(sqrt(3)/4)*area;
end;
	