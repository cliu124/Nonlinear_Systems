function myhist = makehist(tri,x,y,metric)
[edg,edg2tri,tri2edg,nd2tri,nd2edg] = gen_books(tri);
x_ = x(edg); y_ = y(edg);
tx = x_(:,1) -x_(:,2); ty = y_(:,1) -y_(:,2);
c  = [x_(:,1)+x_(:,2) y_(:,1)+y_(:,2)]/2;
t  = [tx  ty ];
tn  = t./repmat(sqrt(t(:,1).^2 +t(:,2).^2 ),1,2);
mm = metric(c);
edgnow  = sum([(tn(:,1) .*mm(:,1) +tn(:,2) .*mm(:,2))  (tn(:,1) .*mm(:,2) +tn(:,2) .*mm(:,3)) ].*tn,2); %ideal edge length at edge center
edgnowl  = sqrt(t(:,1).^2 +t(:,2).^2);
myhist = edgnowl./edgnow;
