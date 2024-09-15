function out = geom_union(in1,in2)
insideA   = and(in1(:,1)>0,in2(:,1)>0);
insideB1  = and(in1(:,1)>=0,in2(:,1)<=0);
insideB2  = and(in1(:,1)<=0,in2(:,1)>=0);
outside   = and(in1(:,1)<0,in2(:,1)<0);
L1 = abs(in1(:,1))<abs(in2(:,1)); L2 = not(L1);
out = zeros(size(in1));
outsideL1 = and(outside,L1);
outsideL2 = and(outside,L2);
out(outsideL1 ,:) = in1(outsideL1,:);
out(outsideL2 ,:) = in2(outsideL2,:);
out(insideB1,:) = in1(insideB1,:);
out(insideB2,:) = in2(insideB2,:);
insideA1 = and(insideA,L1);
insideA2 = and(insideA,L2);
out(insideA1,:) = in2(insideA1,:);
out(insideA2,:) = in1(insideA2,:);