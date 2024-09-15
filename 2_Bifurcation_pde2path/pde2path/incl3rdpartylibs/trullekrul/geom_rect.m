function dist = geom_rect(x_,dxy,poscr)
if nargin == 1
	dxy = repmat(1,1,size(x_,2));
end;
if nargin==2 || nargin == 1
	poscr = repmat(0,1,size(x_,2));
end;
x_ = x_ - repmat(dxy/2+poscr,size(x_,1),1);

LH = dxy(1)/2-abs(x_(:,1));
LV = dxy(2)/2-abs(x_(:,2));
var1 =  x_(:,1)-x_(:,2)-dxy(1)/2+dxy(2)/2;
var2 = -x_(:,1)-x_(:,2)+dxy(1)/2-dxy(2)/2;
var3 =  var1+dxy(1)-dxy(2);
var4 =  var2-dxy(1)+dxy(2);
right = and(var1>0,var2<0); %and(x_(:,1)-dxy(1)/2>x_(:,2)-dxy(2)/2, ...
	   %-x_(:,1)+dxy(1)/2<x_(:,2)+dxy(2)/2);
left  = and(var3<0,var4>0); %and(x_(:,1)+dxy(1)/2<x_(:,2)+dxy(2)/2, ...
	   %-x_(:,1)-dxy(1)/2>x_(:,2)-dxy(2)/2);
uppe = and(0>=var1,var4<=0); %and(x_(:,2)-dxy(2)/2>=x_(:,1)-dxy(1)/2, ...
	  %-x_(:,2)+dxy(2)/2<=x_(:,1)+dxy(1)/2);
lowe = and(0<=var3,var2>=0); %and(x_(:,2)+dxy(2)/2<=x_(:,1)+dxy(1)/2, ...
	  %-x_(:,2)-dxy(2)/2>=x_(:,1)-dxy(1)/2);
inside = and(LV > 0, LH > 0);
chooseH = or(left,right);
chooseV = or(uppe,lowe);
if size(x_,2) == 3
LD = dxy(3)/2-abs(x_(:,3));
var5 =  x_(:,1)-x_(:,3)-dxy(1)/2+dxy(3)/2;
var6 = -x_(:,1)-x_(:,3)+dxy(1)/2-dxy(3)/2;
var7 =  var5+dxy(1)-dxy(3);
var8 =  var6-dxy(1)+dxy(3);
var9  =  x_(:,2)-x_(:,3)-dxy(2)/2+dxy(3)/2;
var10 = -x_(:,2)-x_(:,3)+dxy(2)/2-dxy(3)/2;
var11 =  var9 +dxy(2)-dxy(3);
var12 =  var10-dxy(2)+dxy(3);
right = and(right,and(var5>0,var6<0));
left  = and(left ,and(var7<0,var8>0));
uppe  = and(uppe ,and(var9>0,var10<0));
lowe  = and(lowe ,and(var11<0,var12>0));
front = and(and(var5>=0,var6<=0),and(var9>=0,var10<=0));
behind = and(and(var7<=0,var8>=0),and(var11<=0,var12>=0));
chooseD = or(front,behind);
chooseH = and(or(lowe,uppe),not(chooseD));
chooseV = and(or(right,left),not(chooseD));
%chooseH = or(lowe,uppe);
%chooseV = or(right,left);
chooseD = and(chooseD,not(or(chooseH,chooseV)));
dist = [chooseH.*LH+ ...
        chooseV.*LV+ ...
        chooseD.*LD  ...
        (-sign(x_(:,1)).*chooseH) ... 
        (-sign(x_(:,2)).*chooseV) ...
        (-sign(x_(:,3)).*chooseD)]; %dist vecx vecy vecz
else
dist = [chooseH.*LH+...
        chooseV.*LV ...
        (-sign(x_(:,1)).*chooseH) ... 
        (-sign(x_(:,2)).*chooseV)]; %dist vecx vecy
end;