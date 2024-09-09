function [unconvex,badnds,tri,xy] = mesh_convexity(tri,xy,triQ,bndmesh,Nmetric,options,repeat)
if nargin == 6
	repeat = 10;
end;
for i=1:repeat
unconvex = 0; badnds = [];
nd2tri = inv_table(tri);
if size(xy,2) == 3
	error('not imlpemented');
else %2D
bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
circles = bks_nd2tri2ndO(tri,nd2tri,bndnds_);
NN = sum(circles~=0,2);
[Cc,R] = find(circles'); nbad = size(circles,1);
thmp1 = [0:size(circles,2)-1]'; thmp2 = [2:size(circles,2)+1]';
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
I1 = circles(R+(Cr-1)*nbad); I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
I123 = [I1(:) I2(:) I3(:)];
[bad,z] = elem_inv(I123,xy); %bad(bndnds_(R)) = false; z(bndnds_(R)) = 0;
if any(bad)
	badc = false(size(circles)); badc(R+(Cc-1)*nbad) = bad;
	unconvex = sum(any(badc,2))/size(xy,1);
	[C,I] = min(z);
	badnds = [R(I) I2(I)];
	if options.verbose
		disp(sprintf("%0.0f unconvex out of %0.0f",sum(any(badc,2)),size(xy,1)));
		%figure(); trimesh(tri,xy(:,1),xy(:,2),'color','k');
		%hold on; plot(xy(badnds,1),xy(badnds,2),'ro-'); hold off;
	end;
	delnds_ = false(size(xy,1),1); delnds_(I2(bad)) = true;
        options.fastRM = 0; options.consRM = 0; options.mntn_bks = 1; bks = bks_init(tri);
	bks.rmnd(not(delnds_)) = false;
	[tri,triQ,bks,bndmesh,ndone1,xy,Nmetric] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
else
	break;
end;
end;
end;%fori

close all; [unconvex,badnds,tri_,xy] = mesh_convexity(tri,xy,triQ,bndmesh,Nmetric,options);
tri = sortrows(sort(tri_,2));
trin = sortrows(sort(delaunay(xy(:,1),xy(:,2)),2));
all(tri_(:) == trin(:))
trimesh(tri_,xy(:,1),xy(:,2),'color','k');
hold on; trimesh(trin,xy(:,1),xy(:,2),'color','r'); hold off;