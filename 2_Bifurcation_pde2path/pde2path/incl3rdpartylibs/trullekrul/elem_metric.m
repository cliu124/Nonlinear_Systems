function Xout = elem_metric(tri,xy,p1,p2,p3,p4);
% tri = [1 2 3]; xy = [-0.5 0; 0.5 0; 0 sqrt(3)/2];
max_ratio = 2000; badL = 1e-8; epsdet = 1e-14;%maximum ratio in inverse spac triggers inverse square radii of badL
Nedg = 3;
if nargin == 2 %calculate based
    tris = size(tri,1);
    p1 = xy(tri(:,1),:);
    p2 = xy(tri(:,2),:);
    p3 = xy(tri(:,3),:);
    if size(tri,2) == 4
    	Nedg = 6;
   	p4 = xy(tri(:,4),:);
    end;
else
    tris = size(p1,1);
    if size(p1,2) == 3
    	Nedg = 6;
    end;
end;
r1 = p1-p2; r2 = p1-p3; r3 = p2-p3;
L1 = sqrt(sum(r1.^2,2)); 
L2 = sqrt(sum(r2.^2,2)); 
L3 = sqrt(sum(r3.^2,2));
if nargin == 5 && size(p1,2) == 3 %interior faces!
	Nedg = 3;
	vz = cross(r1,r2,2); 
	vy = cross(vz,r1,2);
	vy = vy./repmat(sqrt(sum(vy.^2,2)),1,3);
	vx = r1./repmat(L1,1,3);
	r1 = [sum(r1.*vx,2) sum(r1.*vy,2)];
	r2 = [sum(r2.*vx,2) sum(r2.*vy,2)];
	r3 = [sum(r3.*vx,2) sum(r3.*vy,2)];
end;
rall = zeros(tris,size(r1,2),Nedg);
if size(r1,2) == 3
	Idm = [ones(tris,1) zeros(tris,1) ones(tris,1) zeros(tris,1) zeros(tris,1) ones(tris,1)];
	r4 = p1-p4; r5 = p2-p4; r6 = p3-p4;
	L4 = sqrt(sum(r4.^2,2)); 
	L5 = sqrt(sum(r5.^2,2)); 
	L6 = sqrt(sum(r6.^2,2)); 
	Lmax = max([L1 L2 L3 L4 L5 L6],[],2); Lmin = min([L1 L2 L3 L4 L5 L6],[],2);
         [Ibad,z] = elem_inv([],[],p1,p2,p3,p4);
         rall(:,:,1) = r1; rall(:,:,2) = r2; rall(:,:,3) = r3;
         rall(:,:,4) = r4; rall(:,:,5) = r5; rall(:,:,6) = r6;
else
	Idm = [ones(tris,1) zeros(tris,1) ones(tris,1)];
	Lmax = max([L1 L2 L3],[],2); Lmin = min([L1 L2 L3],[],2);
	[Ibad,z] = elem_inv([],[],p1,p2,p3);
         rall(:,:,1) = r1; rall(:,:,2) = r2; rall(:,:,3) = r3;
end;
mdet = abs(z);

ratio1 = Lmax.^size(r1,2)./max([mdet epsdet*ones(size(mdet))],[],2);
ratio2 = max([mdet epsdet*ones(size(mdet))],[],2)./Lmin.^size(r1,2);

nIbad = max(ratio1,ratio2)<max_ratio; mdet = mdet(nIbad); fnIbad = find(nIbad);
if any(nIbad)
    %r1 = r1(nIbad,:); r2 = r2(nIbad,:); r3 = r3(nIbad,:); 
    rall = rall(nIbad,:,:);
    tris = nnz(nIbad); %tris - sum(not(nIbad));
    All = zeros(tris,Nedg^2);
    inds = reshape(1:Nedg^2,Nedg,Nedg)';
    for ij=1:Nedg
    	All(:,inds(ij,1)) = rall(:,1,ij).^2; All(:,inds(ij,2)) = 2*rall(:,1,ij).*rall(:,2,ij); All(:,inds(ij,3)) = rall(:,2,ij).^2;
    	if size(r1,2) == 3
    		All(:,inds(ij,4)) = 2*rall(:,1,ij).*rall(:,3,ij); All(:,inds(ij,5)) = 2*rall(:,2,ij).*rall(:,3,ij); All(:,inds(ij,6)) = rall(:,3,ij).^2;
    	end;
    end;
    %     if size(p1,2) == 2
%     A = spalloc(Nedg*tris,Nedg*tris,tris*Nedg^2);
%     C = reshape(permute(reshape(repmat(1:Nedg*tris,Nedg,1),[Nedg,Nedg,tris]),[2,1,3]),[1,Nedg^2*tris]);
%     R = reshape(repmat(1:3*tris,3,1),[1,tris*9]);
%     A(R+(C-1)*3*tris) = reshape([[A11, A12, A13], [A21, A22, A23], [A31, A32, A33]]',[1,9*tris]);
%     end;
    Ain = zeros(Nedg*tris,1+(Nedg-1)*2); 
    ndia = zeros(1,1+(Nedg-1)*2);
    for ij = 1:Nedg
    	for ik=ij:Nedg
    		iks1 = ik:Nedg:Nedg*tris;
    		if ij==1
    		Ain(iks1,ij) = All(:,inds(ik,ik));
    		else
    		iks2 = 1+ik-ij:Nedg:Nedg*tris;
    		Ain(iks1,2*ij-2) = All(:,inds(ik-ij+1,ik));
    		Ain(iks2,2*ij-1) = All(:,inds(ik,ik-ij+1));
    		ndia(2*ij-2) = ij-1;
    		ndia(2*ij-1) = -(ij-1);
		end;
	end;
    end;
    A = spdiags(Ain ,ndia,Nedg*tris,Nedg*tris); %use of spdiags is much faster than linear indexing
   A = A*spdiags(reshape(repmat(1./mdet',Nedg,1),Nedg*tris,1),0,Nedg*tris,Nedg*tris);
    b = ones(tris*Nedg,1);
    X = A\b;
    X = reshape(X',Nedg,tris)'.*repmat(1./mdet,1,Nedg);
    Xout = Idm*badL; Xout(nIbad,:) = X;
else
    tris = size(p1,1); Xout = Idm*badL;
end;

if nargin == 2
	Xout = metric_sqrt(Xout,'invsqrt');
end;
% 
% close all;
% i = 1; N = 101; t = linspace(0,2*pi,N);
% plot(xy(tri(i,:),1)-mean(xy(tri(i,:),1)),xy(tri(i,:),2)-mean(xy(tri(i,:),2)),'-b');
% H = [X(i,1) X(i,2); X(i,2) X(i,3)];
% [v,L] = eig(H); L
% elxy = [cos(t); sin(t)]'*v'*L*v/sqrt(3);
% hold on; plot(elxy(:,1),elxy(:,2),'-r'); hold off;
% areaT = det([r1(i,:); r3(i,:)])/2
% axisprod = prod(1./sqrt(diag(abs(L))))*4
% 
% function X = mysolve(A,b)
% R = 1:3:size(A,1); C = R;
% A11 = full(A(R+0+(C+0-1)*size(A,1)))';
% A12 = full(A(R+0+(C+1-1)*size(A,1)))';
% A13 = full(A(R+0+(C+2-1)*size(A,1)))';
% A21 = full(A(R+1+(C+0-1)*size(A,1)))';
% A22 = full(A(R+1+(C+1-1)*size(A,1)))';
% A23 = full(A(R+1+(C+2-1)*size(A,1)))';
% A31 = full(A(R+2+(C+0-1)*size(A,1)))';
% A32 = full(A(R+2+(C+1-1)*size(A,1)))';
% A33 = full(A(R+2+(C+2-1)*size(A,1)))';


% tris = 4;
% A11 =   ones(tris,1); A12 = 2*ones(tris,1); A13 = 3*ones(tris,1);
% A21 = 4*ones(tris,1); A22 = 5*ones(tris,1); A23 = 6*ones(tris,1);
% A31 = 7*ones(tris,1); A32 = 8*ones(tris,1); A33 = 9*ones(tris,1);
%  [[A11, A12, A13], [A21, A22, A23], [A31, A32, A33]]
 
%[xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(20,gen_options(),0);  Xout = elem_metric(tri,xy); [eigL,eigR] = analyt_eig(Xout); [tmp,vol] = elem_inv(tri,xy); theE = abs(1-prod(eigL,2)*sqrt(2)/12./(vol/6)); sum(theE>1e-6)
%.^-0.5
%i=4; [t1,t2] = eig([Xout(i,1) Xout(i,2) Xout(i,4); Xout(i,2) %Xout(i,3) Xout(i,5); Xout(i,4) Xout(i,5) Xout(i,6)]);
%t2 = diag(diag(t2).^-0.5);
%[abs(1-prod(diag(t2))*sqrt(2)/12./(vol(i)/6)) theE(i)]

%[xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(20,gen_options(),1);  Xout = elem_metric(tri,xy); [eigL,eigR] = analyt_eig(Xout); [tmp,vol] = elem_inv(tri,xy); theE = abs(1-prod(eigL,2).^-0.5*sqrt(2)/12./(vol/6)); sum(theE>1e-6)