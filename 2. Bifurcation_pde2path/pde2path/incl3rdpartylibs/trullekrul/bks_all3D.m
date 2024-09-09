function [fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac,nd2edg,nd2tri] = bks_all3D(tri);
%gen faces
faca = [tri(:,1) tri(:,2) tri(:,4); ...
        tri(:,2) tri(:,1) tri(:,3); ...
        tri(:,3) tri(:,4) tri(:,2); ...
        tri(:,4) tri(:,3) tri(:,1)];
faca2tri = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1) 1:size(tri,1)]';

[fac,I] = sort(faca,2); %biggest in last column, smaller in first
I = even_permut(I);
[fac,Is] = sortrows(fac);
tris = faca2tri(Is);
I = I(Is);
d =[true; or(or(fac(1:end-1,1)~=fac(2:end,1), ...
	       fac(1:end-1,2)~=fac(2:end,2)), ...
	       fac(1:end-1,3)~=fac(2:end,3))]; Nd = sum(d);
ad = not(d); Nad = sum(ad);
d2 = [ad(2:end); false]; d2 = d2(d);
fac = fac(d,:);

% gen fac2tri
tris2 = zeros(Nd,1); tris2(d2) = tris(ad); 
fac2tri = [tris(d) tris2];

%gen fac2tri2, the element to which the permutation of the face belongs 
fac2tri2 = faca2tri(Is(d));
%fac2tri2 = fac2tri(:,1); fac2tri2(i) = fac2tri(i,2);
%flip fac2tri according to fac permutation
i = and(d2,I(d)); fac2tri(i,:) = fac2tri(i,[2 1]); 

%gen faca2fac
faca2fac = zeros(size(faca,1),1);
faca2fac(Is(d))  = 1:size(fac,1);
faca2fac(Is(ad)) = find(d2);
% gen tri2faca, tri2fac
tri2faca = reshape([1:size(tri,1)*4]',size(tri,1),4);
tri2fac = faca2fac(tri2faca);

if nargout <= 3
	return;
end;

%gen edges
edga = [faca(:,2) faca(:,1); ...
        faca(:,3) faca(:,2); ...
        faca(:,1) faca(:,3)];
edga2faca = [1:size(faca,1) 1:size(faca,1) 1:size(faca,1)]';
edga2tri  = [faca2tri; faca2tri; faca2tri];
faca2edga = reshape(1:size(edga,1),size(faca,1),3);
fac2edga = faca2edga(Is(d),:);

edg = sort(edga,2); %biggest in last column
[edg,Ise] = sortrows(edg);
facs = edga2faca(Ise);
tris = edga2tri(Ise);
de =[true; or(edg(1:end-1,1)~=edg(2:end,1),...
              edg(1:end-1,2)~=edg(2:end,2))]; Nde = sum(de);
ade = not(de); Nad = sum(ade);
de2 = cumsum(de); %de2 = [ade(2:end); false]; de2 = de2(de);

edg = edg(de,:);
%gen edg2faca
edgs = cumsum(de); %edgs = zeros(size(edga,1),1); edgs(de) = 1; edgs = cumsum(edgs);
edg2faca = rpval2M(edgs,facs); edg2faca = rpval2M_clean(edg2faca);
%gen edg2fac
edg2faca_ = edg2faca; edg2faca_(edg2faca_==0) = numel(faca2fac)+1;
faca2fac_ = [faca2fac; 0];
edg2fac = faca2fac_(edg2faca_); edg2fac = rpval2M_clean(edg2fac);
%gen edg2tri
edg2tri = rpval2M(edgs,tris); edg2tri = rpval2M_clean(edg2tri);

%gen edga2edg
edga2edg = zeros(size(edga,1),1);
edga2edg(Ise(de))  = 1:size(edg,1);
edga2edg(Ise(ade)) = de2(ade); %find(de2);
% gen tri2edga, tri2edg
tri2edga = reshape([1:size(tri,1)*12]',size(tri,1),12);
tri2edg = reshape(edga2edg(tri2edga(repmat([true false true false true true false false true false true false],size(tri,1),1))),size(tri,1),6);
%tri2edg = sort(edga2edg(tri2edga)');
%tri2edg = reshape(tri2edg([true(1,size(tri,1); tri2edg(1:end-1,:)~= tri2edg(2:end,:)]),6,size(tri,1))';

%gen faca2edg, fac2edg (fac2edga does not make sense)
faca2edg = edga2edg(faca2edga);
fac2edg = faca2edg(Is(d),:);

if nargout == 21
	return;
end;

% gen nd2fac
nd2fac = inv_table(fac);
%fact = fac;
%[nds,I] = sort(fact(:));
%facs = [1:size(fac,1) 1:size(fac,1) 1:size(fac,1)]'; facs = facs(I);
%nd2fac =  rpval2M(nds,facs);

if nargout == 22
	return;
end;
% gen nd2edg
nd2edg = inv_table(edg);
%edgt = edg;
%[nds,I] = sort(edgt(:));
%edgs = [1:size(edg,1) 1:size(edg,1)]'; edgs = edgs(I);
%nd2edg =  rpval2M(nds,edgs);

if nargout == 23
	return;
end;
% gen nd2tri
nd2tri = inv_table(tri);
%trit = tri;
%[nds,I] = sort(trit(:));
%tris = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1)  1:size(tri,1)]'; tris = tris(I);
%nd2tri = rpval2M(nds,tris);

function is_even = even_permut(I)
is_even = false(size(I,1),1);
is_even(all(I == repmat([1,2,3],size(I,1),1),2)) = true;
is_even(all(I == repmat([3,1,2],size(I,1),1),2)) = true;
is_even(all(I == repmat([2,3,1],size(I,1),1),2)) = true;
%[1 3 2]: 0
%[1 2 3]: 1
%[2 1 3]: 0
%[2 3 1]: 1
%[3 2 1]: 0
%[3 1 2]: 1

