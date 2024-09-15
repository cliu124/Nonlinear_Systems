function Nmetric = metric_pnorm(tri,xy,z,eta,pnorm,nd2tri)
%p-norm scaling to the metric, as in Chen, Sun and Xu, Mathematics of
%  Computation, Volume 76, Number 257, January 2007, pp. 179-204.

%z = sin(xy(:,1)).*cos(xy(:,2)); trixy = [mean(reshape(xy(tri,1),size(tri)),2) mean(reshape(xy(tri,2),size(tri)),2)]; gradz = [cos(trixy(:,1)).*cos(trixy(:,2)) -sin(trixy(:,1)).*sin(trixy(:,2))]; gradzn = [cos(xy(:,1)).*cos(xy(:,2)) -sin(xy(:,1)).*sin(xy(:,2))]; heszn = [-sin(xy(:,1)).*cos(xy(:,2)) -cos(xy(:,1)).*sin(xy(:,2)) -cos(xy(:,1)).*sin(xy(:,2)) -sin(xy(:,1)).*cos(xy(:,2))]; hesz = [-sin(trixy(:,1)).*cos(trixy(:,2)) -cos(trixy(:,1)).*sin(trixy(:,2)) -cos(trixy(:,1)).*sin(trixy(:,2)) -sin(trixy(:,1)).*cos(trixy(:,2))];
if nargin == 5 && numel(pnorm) == 1
	nd2tri = inv_table(tri);
elseif nargin == 5 
    nd2tri = pnorm; pnorm = 2;
elseif nargin == 4
    pnorm = 2; nd2tri = inv_table(tri);
end;
nd2tri_ = nd2tri~=0;
dim = size(xy,2);
xys = size(xy,1);
tris = size(tri,1);
triz = z(tri(:,2:end)) - repmat(z(tri(:,1)),1,dim);
b = reshape(permute(reshape(repmat([1; repmat([zeros(dim,1); 1],dim-1,1)],1,tris),[dim dim tris]),[1 3 2]),tris*dim,dim);
%R = reshape(repmat(1:size(xy,1),sqdim,1),size(b));
%C = repmat(1:size(xy,2)
if size(xy,2) == 2
	det11 = xy(tri(:,2),1)-xy(tri(:,1),1);
	det12 = xy(tri(:,3),1)-xy(tri(:,1),1);
	det21 = xy(tri(:,2),2)-xy(tri(:,1),2);
	det22 = xy(tri(:,3),2)-xy(tri(:,1),2);
	All = [det11 det12 det21 det22]; dem = det11.*det22-det21.*det12;
else %3D
	det11 = xy(tri(:,2),1)-xy(tri(:,1),1);
	det12 = xy(tri(:,3),1)-xy(tri(:,1),1);
	det13 = xy(tri(:,4),1)-xy(tri(:,1),1);
	det21 = xy(tri(:,2),2)-xy(tri(:,1),2);
	det22 = xy(tri(:,3),2)-xy(tri(:,1),2);
	det23 = xy(tri(:,4),2)-xy(tri(:,1),2);
	det31 = xy(tri(:,2),3)-xy(tri(:,1),3);
	det32 = xy(tri(:,3),3)-xy(tri(:,1),3);
	det33 = xy(tri(:,4),3)-xy(tri(:,1),3);
	All = [det11 det12 det13 det21 det22 det23 det31 det32 det33];
	dem = det11.*det22.*det33 + det12.*det23.*det31 + det13.*det21.*det32 ...
	     -det21.*det12.*det33 - det11.*det23.*det32 - det22.*det13.*det31;
end;
nd2trim = zeros(size(nd2tri)); 
nd2trim(nd2tri_) = dem(nd2tri(nd2tri_));
nd2trim = nd2trim./repmat(sum(nd2trim,2),1,size(nd2tri,2));
%Ain = zeros(size(tri,1)*dim,dim*2-1);
%ndia = zeros(1,size(Ain,2));
%inds = reshape(1:dim^2,dim,dim)';
%for ij = 1:dim
    	%for ik=ij:dim
    		%iks1 = ik:dim:dim*tris;
    		%if ij==1
    		%Ain(iks1,ij) = All(:,inds(ik,ik));
    		%else
    		%iks2 = 1+ik-ij:dim:dim*tris;
    		%Ain(iks1,2*ij-2) = All(:,inds(ik-ij+1,ik));
    		%Ain(iks2,2*ij-1) = All(:,inds(ik,ik-ij+1));
    		%ndia(2*ij-2) = ij-1;
    		%ndia(2*ij-1) = -(ij-1);
		%end;
	%end;
%end;
%tic;
%A = spdiags(Ain,ndia,size(b,1),size(b,1));
%toc
%[R,C,vals] = find(A);
All = reshape(All',numel(All),1);
R = reshape(repmat(1:tris*dim,dim,1),numel(All),1);
C = reshape(permute(reshape(R,[dim dim tris]),[2 1 3]),numel(All),1); 
A = sparse(R,C,All,tris*dim,tris*dim);
X = A\b;
X = reshape(X',dim^2,tris)';
gradf = zeros(size(triz));
%gradfn = zeros(size(xy,1),dim);
hes = zeros(tris,dim*dim);
hesn = zeros(xys,dim*dim);
for i=1:dim
for j=1:dim
	gradf(:,i) = gradf(:,i) + triz(:,j).*X(:,dim*(j-1)+i);
end;
gradfn_ = zeros(size(nd2tri)); gradfn_(nd2tri_) = gradf(nd2tri(nd2tri_),i);
%gradfn(:,i) = sum(gradfn_.*nd2trim,2);
gradfn_ = sum(gradfn_.*nd2trim,2);
triz_ = gradfn_(tri(:,2:end)) - repmat(gradfn_(tri(:,1)),1,dim);
for j=1:dim
for k=1:dim
	hes(:,(i-1)*dim+j) = hes(:,(i-1)*dim+j) + triz_(:,k).*X(:,dim*(k-1)+j);
end;
hesn_ = zeros(size(nd2tri)); hesn_(nd2tri_) = hes(nd2tri(nd2tri_),(i-1)*dim+j);
hesn(:,(i-1)*dim+j) = sum(hesn_.*nd2trim,2);
end;
end;

if size(xy,2) == 2
	hesn(:,2) = mean(hesn(:,2:3),2);
	hesn = hesn(:,[1 2 4]);
    hesn(:,[1 3]) = hesn(:,[1 3]) + 1e-5;
    %detn = hesn(:,1).*hesn(:,3)-hesn(:,2).^2.;
else
	hesn(:,2) = mean(hesn(:,[2 4]),2);
	hesn(:,3) = mean(hesn(:,[3 7]),2);
	hesn(:,5) = mean(hesn(:,[5 8]),2);
	hesn = hesn(:,[1 2 5 3 5 9]);
    hesn(:,[1 3 6]) = hesn(:,[1 3 6]) + 1e-5;
    %detn = hesn(:,1).*hesn(:,3).*hesn(:,6)+2.*hesn(:,2).*hesn(:,4).*hesn(:,5) ...
    %     - hesn(:,6).*hesn(:,2).^2. - hesn(:,1).*hesn(:,5).^2. - hesn(:,3).*hesn(:,4).^2.;
end;
[eigL,eigR] = analyt_eig(hesn);
eigL = reshape(max(abs(eigL),1e-12),size(eigL));
detn = prod(eigL,2);
%hesn = analyt_prod(analyt_fulleig(eigL),eigR);
hesn = analyt_prod(analyt_fulleig(eigL.^0.5),eigR);
%exponent = -1./(2.*pnorm + size(xy,2));
exponent = -0.5/(2.*pnorm + size(xy,2));
%Nmetric = 1/eta*repmat(sqrt(detn).^exponent,1,size(hesn,2)).*hesn;
Nmetric = 1/sqrt(eta)*repmat(detn.^exponent,1,size(hesn,2)).*hesn;