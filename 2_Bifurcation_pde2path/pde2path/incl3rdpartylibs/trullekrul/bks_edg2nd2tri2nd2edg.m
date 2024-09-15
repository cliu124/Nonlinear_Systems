function [nghOa,ngh,deltri] = bks_edg2nd2tri2nd2edg(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,options,nd2fac)
tris = size(tri,1);
edgs = size(edg,1);
edgas = max(edga2edg);
edgas2 = size(edga2edg,1);
nds = size(nd2tri,1);

%compute deltri
nd2tri_ = nd2tri; nd2tri_(nd2tri_==0) = tris + 1;
deltri = reshape(nd2tri_(edg(:),:),edgs,size(nd2tri,2)*2); 
deltri = sort(deltri,2); I = and(deltri~=tris+1,[true(edgs,1) deltri(:,2:end)~=deltri(:,1:end-1)]); deltri(not(I)) = tris+1;
deltri_ = sort(deltri,2); maxdeltri = max(sum(deltri_~=tris+1,2));
deltri_=deltri_(:,1:maxdeltri); deltri=deltri_; deltri(deltri==tris+1) = 0;
%generate ngh
tri_ = [tri; repmat(tris+1,1,size(tri,2))];
trinds = reshape(tri_(deltri_(:),:),edgs,size(tri_,2)*size(deltri_,2));
[trinds,tmp] = sort(trinds,2);
trinds([false(edgs,1) trinds(:,2:end)==trinds(:,1:end-1)]) = nds+1;
[trinds,tmp] = sort(trinds,2);
maxtrinds = max(sum(trinds~=nds+1,2)); trinds = trinds(:,1:maxtrinds);
nd2edg_ = [nd2edg; repmat(0,1,size(nd2edg,2))]; nd2edg_(nd2edg_==0) = edgs+1;
ngh = reshape(nd2edg_(trinds(:),:),edgs,size(nd2edg_,2)*size(trinds,2));
[ngh,tmp] = sort(ngh,2);
I = or([false(edgs,1) ngh(:,2:end)==ngh(:,1:end-1)],repmat((1:edgs)',1,size(ngh,2))==ngh);
ngh(I) = edgs+1;
[ngh,tmp] = sort(ngh,2);
maxngh = max(sum(ngh~=edgs+1,2)); 
ngh = ngh(:,1:maxngh);
ngh(ngh==edgs+1) = 0;

%generate nghOa
tri2edga_ = [tri2edga; repmat(edgas2+1,1,size(tri2edga,2))];
nghsa = reshape(tri2edga_(deltri_(:),:),edgs,size(deltri_,2)*size(tri2edga_,2));
edga2edg_ = [edga2edg; edgas+1];
[nghs,C] = sort(edga2edg_(nghsa)');
In = or([false(1,edgs); nghs(1:end-1,:)==nghs(2:end,:)],[nghs(1:end-1,:)==nghs(2:end,:); false(1,edgs)]); 
%taking out doubles
I = and(nghs~=edgas+1,not(In)); R = repmat(1:edgs,size(nghs,1),1);
Ina = false(size(nghsa)); Ina(R(In)+(C(In)-1)*size(Ina,1)) = true;
if nargin == 8
ntrnl = reshape(nd2edg_(edg(:),:),edgs,size(nd2edg_,2)*2)'; inds = [2 1]; %[1 2];
else
nd2fac_ = nd2fac; nd2fac_(nd2fac_==0) = edgas+1;
ntrnl = reshape(nd2fac_(edg(:),:),edgs,size(nd2fac_,2)*2)'; inds = [2 1];
end; 
if options.nosparse
tbl1 = [nghs(I) R(I) ]; nvtbl1 = inv_table(tbl1(:,inds(1))); delC = C(I);
[ntrnl,tmp] = sort(ntrnl);
I = and(ntrnl~=edgas+1,[true(1,edgs); ntrnl(1:end-1,:)~=ntrnl(2:end,:)]);
R = repmat(1:edgs,size(ntrnl,1),1);
tbl2 = [ntrnl(I) R(I)];
tbl2tbl =  bks_bndedg2edg(tbl1(:,inds),nvtbl1,tbl2(:,inds));
R = tbl1(tbl2tbl,2); delC = delC(tbl2tbl);
else
A = sparse(nghs(I),R(I),C(I),edgas,edgs);
[ntrnl,tmp] = sort(ntrnl);
I = and(ntrnl~=edgas+1,[true(1,edgs); ntrnl(1:end-1,:)~=ntrnl(2:end,:)]);
R = repmat(1:edgs,size(ntrnl,1),1);
A2 = sparse(ntrnl(I),R(I),true(nnz(I(:)),1),edgas,edgs);
A = A.*A2; [tmp,R,delC] = find(A); clear A2 A;%it is this very product that is ''tricky'' to do using sorts
end;
Ina(R+(delC-1)*size(Ina,1)) = true;
%R = repmat(1:edgs,size(nghsa,2),1);
%nghsa(R(In)+(C(In)-1)*size(nghsa,1)) = edgas2+1;
nghsa(Ina) = edgas2+1;
nghOa = sort(nghsa,2);
maxnghs = max(sum(nghOa~=edgas2+1,2)); nghOa = nghOa(:,1:maxnghs);
nghOa(nghOa==edgas2+1) = 0;

%function [nghOa,ngh,deltri] = ngh_edg2(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,nd2fac)
%if nargin == 8
%chnksz = 5000;
%else %2D
%chnksz = 50000;
%end;
%bmin = 1:chnksz:size(edg,1);
%bmax = [chnksz:chnksz:size(edg,1) size(edg,1)];
%nghOa = []; ngh = []; deltri = [];
%for i=1:numel(bmin)
	%edgi = bmin(i):bmax(i);
	%if nargin == 8
	%[nghOa_,ngh_,deltri_] = ngh_edg2a(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,edgi,nd2fac);
	%else
	%[nghOa_,ngh_,deltri_] = ngh_edg2a(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,edgi);
	%end;
%deltri = [deltri zeros(size(deltri,1),size(deltri_,2)-size(deltri,2))]; 
%nghOa  = [nghOa  zeros(size(nghOa ,1),size(nghOa_ ,2)-size(nghOa ,2))]; 
%ngh    = [ngh    zeros(size(ngh   ,1),size(ngh_   ,2)-size(ngh   ,2))]; 
%deltri_ = [deltri_ zeros(size(deltri_,1),size(deltri,2)-size(deltri_,2))]; 
%nghOa_  = [nghOa_  zeros(size(nghOa_ ,1),size(nghOa ,2)-size(nghOa_ ,2))]; 
%ngh_    = [ngh_    zeros(size(ngh_   ,1),size(ngh   ,2)-size(ngh_   ,2))]; 
%deltri = [deltri; deltri_];
%nghOa  = [nghOa ; nghOa_ ];
%ngh    = [ngh   ; ngh_   ];
%end;
%
%function [nghOa,ngh,deltri] = ngh_edg2a(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,edgi,nd2fac)
%tris = size(tri,1);
%edgs = numel(edgi);
%edgs2 = size(edg,1);
%edgas = max(edga2edg);
%edgas2 = size(edga2edg,1);
%nds = size(nd2tri,1);
%
%%compute deltri
%nd2tri_ = nd2tri; nd2tri_(nd2tri_==0) = tris + 1;
%deltri = reshape(nd2tri_(edg(edgi,:),:),edgs,size(nd2tri,2)*2); 
%deltri = sort(deltri,2); I = and(deltri~=tris+1,[true(edgs,1) deltri(:,2:end)~=deltri(:,1:end-1)]); deltri(not(I)) = tris+1;
%deltri_ = sort(deltri,2); maxdeltri = max(sum(deltri_~=tris+1,2));
%deltri_=deltri_(:,1:maxdeltri); deltri=deltri_; deltri(deltri==tris+1) = 0;
%%generate ngh
%tri_ = [tri; repmat(tris+1,1,size(tri,2))];
%trinds = reshape(tri_(deltri_(:),:),edgs,size(tri_,2)*size(deltri_,2));
%[trinds,tmp] = sort(trinds,2);
%trinds([false(edgs,1) trinds(:,2:end)==trinds(:,1:end-1)]) = nds+1;
%[trinds,tmp] = sort(trinds,2);
%maxtrinds = max(sum(trinds~=nds+1,2)); trinds = trinds(:,1:maxtrinds);
%nd2edg_ = [nd2edg; repmat(0,1,size(nd2edg,2))]; nd2edg_(nd2edg_==0) = edgs2+1;
%ngh = reshape(nd2edg_(trinds(:),:),edgs,size(nd2edg_,2)*size(trinds,2));
%[ngh,tmp] = sort(ngh,2);
%I = or([false(edgs,1) ngh(:,2:end)==ngh(:,1:end-1)],repmat(edgi',1,size(ngh,2))==ngh);
%ngh(I) = edgs2+1;
%[ngh,tmp] = sort(ngh,2);
%maxngh = max(sum(ngh~=edgs2+1,2)); 
%ngh = ngh(:,1:maxngh);
%ngh(ngh==edgs2+1) = 0;
%
%%generate nghOa
%tri2edga_ = [tri2edga; repmat(edgas2+1,1,size(tri2edga,2))];
%nghsa = reshape(tri2edga_(deltri_(:),:),edgs,size(deltri_,2)*size(tri2edga_,2));
%edga2edg_ = [edga2edg; edgas+1];
%[nghs,C] = sort(edga2edg_(nghsa)');
%In = or([false(1,edgs); nghs(1:end-1,:)==nghs(2:end,:)],[nghs(1:end-1,:)==nghs(2:end,:); false(1,edgs)]); %taking out doubles
%I = and(nghs~=edgas+1,not(In)); R = repmat(1:edgs,size(nghs,1),1);
%Ina = false(size(nghsa)); Ina(R(In)+(C(In)-1)*size(Ina,1)) = true;
%A = sparse(nghs(I),R(I),C(I),edgas,edgs);
%if nargin == 8
%ntrnl = reshape(nd2edg_(edg(edgi,:),:),edgs,size(nd2edg_,2)*2)';
%else
%nd2fac_ = nd2fac; nd2fac_(nd2fac_==0) = edgas+1;
%ntrnl = reshape(nd2fac_(edg(edgi,:),:),edgs,size(nd2fac_,2)*2)';
%end;
%[ntrnl,tmp] = sort(ntrnl);
%I = and(ntrnl~=edgas+1,[true(1,edgs); ntrnl(1:end-1,:)~=ntrnl(2:end,:)]);
%R = repmat(1:edgs,size(ntrnl,1),1);
%A2 = sparse(ntrnl(I),R(I),true(nnz(I(:)),1),edgas,edgs);
%A = A.*A2; [tmp,R,delC] = find(A); clear A2 A;%it is this very product that is ''tricky'' to do using sorts
%Ina(R+(delC-1)*size(Ina,1)) = true;
%%R = repmat(1:edgs,size(nghsa,2),1);
%%nghsa(R(In)+(C(In)-1)*size(nghsa,1)) = edgas2+1;
%nghsa(Ina) = edgas2+1;
%nghOa = sort(nghsa,2);
%maxnghs = max(sum(nghOa~=edgas2+1,2)); nghOa = nghOa(:,1:maxnghs);
%nghOa(nghOa==edgas2+1) = 0;



%function [nghOa,ngh,deltri] = ngh_edg2(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,options);
%tris = size(tri,1);
%edgs = size(edg,1);
%nds = size(nd2tri,1);
%
%nd2tri_ = nd2tri; nd2tri_(nd2tri_==0) = tris + 1;
%tri2edga_ = [tri2edga; [0 0 0]];
%tri2edga_1 = tri2edga_(:,1); tri2edga_2 = tri2edga_(:,2); tri2edga_3 = tri2edga_(:,3);
%deltri = [nd2tri_(edg(:,1),:) nd2tri_(edg(:,2),:)];
%deltri = sort(deltri,2); I = and(deltri~=tris+1,[true(edgs,1) deltri(:,2:end)~=deltri(:,1:end-1)]); deltri(not(I)) = tris+1;
%deltri_ = sort(deltri,2); deltri=deltri_; deltri(deltri==tris+1) = 0;
%nghsa = [tri2edga_1(deltri_) tri2edga_2(deltri_) tri2edga_3(deltri_)]; clear tri2edga_1 tri2edga_2 tri2edga_3;
%nghsa_ = nghsa; nghsa_(nghsa_==0) = size(edga,1)+1; edga2edg_ = [edga2edg; 0];
%%nghs = edga2edg_(nghsa_);
%%we actually need edg->nd->tri->nd->edg
%R = repmat((1:edgs),size(deltri_,2),1);
 %deltri_ = deltri_'; R = R(deltri_~=tris+1); deltri_ = deltri_(deltri_~=tris+1);
%tri1 = tri(deltri_,1); tri2 = tri(deltri_,2); tri3 = tri(deltri_,3);
%nghedg = [nd2edg(tri1,:)'; nd2edg(tri2,:)'; nd2edg(tri3,:)']; clear tri1 tri2 tri3;
%R = repmat(R',size(nghedg,1),1); I = nghedg~=0; R = R(I); nghedg = nghedg(I);
%nghs = rpval2M(R,nghedg);
%
%% take out doubles and self, gen ngh
%nghsS = sort(nghs,2); 
%I = and(nghsS~=repmat((1:edgs)',1,size(nghs,2)),and(nghsS~=0,[nghsS(:,1:end-1) ~=nghsS(:,2:end) true(edgs,1)]));
%[C,R] = find(I'); nI = not(I); C_ = nghsS(R+(C-1)*size(nghs,1));
%ngh = rpval2M(R,C_);
%
%%take out interior edges as doubles as well as using nd2edg, gen nghOa
%[nghsS,Co] = sort(edga2edg_(nghsa_),2); Ro = repmat((1:edgs)',1,size(nghsS,2));
%I = and(nghsS~=repmat((1:edgs)',1,size(nghsS,2)),and(nghsS~=0,and([nghsS(:,1:end-1) ~=nghsS(:,2:end) true(edgs,1)],[true(edgs,1) nghsS(:,1:end-1) ~=nghsS(:,2:end)])));
%[C,R] = find(I'); nI = not(I);
%nghsa(Ro(nI)+(Co(nI)-1)*edgs) = 0;
%
%nghsa = sort(nghsa,2); columns = max(sum(nghsa~=0,2)); nghsa = nghsa(:,size(nghsa,2):-1:size(nghsa,2)-columns+1);
%nghsa_ = nghsa; nghsa_(nghsa_==0) = size(edga,1)+1; nghs = edga2edg_(nghsa_);
%[nghsS,Co] = sort(nghs,2); Ro = repmat((1:edgs),size(nghs,2),1); nghsS = nghsS'; Co = Co';
%I = nghsS~=0; nI = not(I);
%[C,R] = find(I); 
%A = sparse(nghsS(I),R,C);
%nd2edga = [nd2edg(edg(:,1),:) nd2edg(edg(:,2),:)]; 
%R_ = repmat((1:edgs),size(nd2edga,2),1); nd2edga = nd2edga'; nd2edga(nd2edga==R_) = 0;
%I_= nd2edga ~= 0;
%A2 = sparse(nd2edga(I_),R_(I_),ones(sum(I_(:)),1));
%Am = A2.*A; [Cm,Rm,Im] = find(Am); clear Am A2 A;%it is this very product that is ''tricky'' to do using sorts
%nI(Im+(Rm-1)*size(nI,1)) = true;
%nghsa(Ro(nI)+(Co(nI)-1)*edgs) = 0;
%nghOa = sort(nghsa,2); columns = max(sum(nghOa~=0,2)); nghOa = nghOa(:,size(nghOa,2):-1:size(nghOa,2)-columns+1);