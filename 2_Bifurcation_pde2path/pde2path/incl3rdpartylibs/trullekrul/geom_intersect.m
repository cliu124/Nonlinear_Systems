function out = geom_intersect(in1,in2)
outsideA   = and(in1(:,1)<0,in2(:,1)<0);
outsideB1  = and(in1(:,1)>0,in2(:,1)<0);
outsideB2  = and(in1(:,1)<0,in2(:,1)>0);
inside   = and(in1(:,1)>=0,in2(:,1)>=0);
L1 = abs(in1(:,1))<abs(in2(:,1)); L2 = not(L1);
out = zeros(size(in1));
insideL1 = and(inside,L1);
insideL2 = and(inside,L2);
out(insideL1 ,:) = in1(insideL1,:);
out(insideL2 ,:) = in2(insideL2,:);
out(outsideB1,:) = in2(outsideB1,:);
out(outsideB2,:) = in1(outsideB2,:);
outsideA1 = and(outsideA,L1);
outsideA2 = and(outsideA,L2);
out(outsideA1,:) = in1(outsideA1,:);
out(outsideA2,:) = in2(outsideA2,:);