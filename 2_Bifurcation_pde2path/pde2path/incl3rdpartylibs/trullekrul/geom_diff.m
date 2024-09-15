function out = geom_diff(in1,in2)
inside   = and(in1(:,1)>=0,in2(:,1)<=0);
outsideA = and(in1(:,1)< 0,in2(:,1) <0);
outsideB = in2(:,1)>0;
L1 = abs(in1(:,1))<abs(in2(:,1)); L2 = not(L1);
out = zeros(size(in1));
insideL1 = and(inside,L1);
insideL2 = and(inside,L2);
out(insideL1 ,:) =  in1(insideL1,:);
out(insideL2 ,:) = -in2(insideL2,:);
out(outsideA,:) =  in1(outsideA,:);
out(outsideB,:) = -in2(outsideB,:);
%insideA1 = and(insideA,L1);
%insideA2 = and(insideA,L2);
%out(insideA1,:) = in2(insideA1,:);
%out(insideA2,:) = in1(insideA2,:);