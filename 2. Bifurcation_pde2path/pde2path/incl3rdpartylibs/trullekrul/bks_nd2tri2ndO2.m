% ALTERNATIVE IN TERMS OF (ordered) EDGE NUMBERS
function [ngh,nghO] = bks_nd2tri2ndO2(edg,nd2edg,nd2tri,edg2tri,edga,tri2edga,edga2edg);
tris = size(tri2edga,1);
edgs = size(edg2tri,1);
nds = size(nd2tri,1);
bndedg = sum(edg2tri>0,2)==1;
bndnds_ = sort(reshape(edg(bndedg,:),sum(bndedg)*2,1));
bndnds = false(nds,1); bndnds(bndnds_) = true; bndnds_ = [bndnds; false];
nbndnds = true(nds,1); nbndnds(bndnds) = false;
maxnghI  = max(sum(nd2tri(nbndnds,:)~=0,2));
maxnghb = max(sum(nd2tri(bndnds,:)~=0,2))+1;
maxngh = max(maxnghI,maxnghb);

nd2tri_ = nd2tri; nd2tri_(nd2tri_==0) = tris + 1;
tri2edga_ = [tri2edga; [0 0 0]];
tri2edga_1 = tri2edga_(:,1); tri2edga_2 = tri2edga_(:,2); tri2edga_3 = tri2edga_(:,3);
nghsa = [tri2edga_1(nd2tri_) tri2edga_2(nd2tri_) tri2edga_3(nd2tri_)];
nghsa_ = nghsa; nghsa_(nghsa_==0) = size(edga,1)+1; edga2edg_ = [edga2edg; 0];
nghs = edga2edg_(nghsa_);
%%delete boundary edges (for boundary nodes)
%bndedg_ = [(1:edgs)'; 0]; bndedg_(bndedg) = 0;
%nghs_ = nghs(bndnds,:); nghs_(nghs_==0) = edgs+1;
%nghs(bndnds,:)  = bndedg_(nghs_);
%%delete "interior edges"
%[nghs, C] = sort(nghs,2); R=repmat((1:nds)',1,size(nghs,2));
%nghsa = nghsa(R+(C-1)*nds);
%IO = and(and([true(nds,1) nghs(:,2:end) ~= nghs(:,1:end-1)],[nghs(:,2:end) ~= nghs(:,1:end-1) true(nds,1)]),nghs~=0);
%nghO = zeros(size(nghs)); nghO(IO) = nghsa(IO);
%%clean out zeros
%columns = max(sum(nghO~=0,2));
%nghO = sort(nghO,2); nghO = nghO(1:nds,size(nghO,2)-columns+1:end);

%delete interior edges using nd2edg and sparse()
%R1 = repmat((1:nds)',1,size(nghs,2)); I = nghs~=0; C1=nghs(I); R1=R1(I);
%A1 = sparse(R1,C1,ones(size(R1))); A1a = sparse(R1,C1,nghsa(I));
%R2 = repmat((1:nds)',1,size(nd2edg,2)); I = nd2edg~=0; C2=nd2edg(I); R2=R2(I);
%A2 = sparse(R2,C2,ones(size(R2)));
%A = A1+A2; [C_,R] = find(A'==1); A=sparse(C_,R,ones(size(R))); A=A1a'.*A; [C,R,C_] = find(A);

R1 = repmat((1:nds),size(nghs,2),1); nghst=nghs'; nghsat = nghsa'; I = nghst~=0; C1=nghst(I); R1=R1(I);
A1 = sparse(C1,R1,ones(size(R1))); A1a = sparse(C1,R1,nghsat(I));
R2 = repmat((1:nds),size(nd2edg,2),1); nd2edgt = nd2edg'; I = nd2edgt~=0; C2=nd2edgt(I); R2=R2(I);
A2 = sparse(C2,R2,ones(size(R2)));
A = A1+A2; [C_,R] = find(A==1); A=sparse(C_,R,ones(size(R))); A=A1a.*A; [C,R,C_] = find(A);

nN = diff([0 find([R(1:end-1)'~=R(2:end)' true])]); columns = max(nN);
C = repmat((1:max(nN))',1,numel(nN)); CnN = repmat(nN,max(nN),1); I = C <= CnN; 
nghO = zeros(numel(nN),columns); nghO(R+(C(I)-1)*numel(nN)) = C_;
nghO = nghO(:,size(nghO,2):-1:1);

%ngh = edga2circ(nghO);
%function ngh = edga2circ(nghO,edga);
edga_ = [edga; [0 0]];
nghO_ = nghO; nghO_(nghO_==0) = size(edga_,1);
nghO_ = edga_(reshape(nghO_',numel(nghO_),1),:);
%set first nodes
ngh = zeros(nds,maxngh);
ngh(nbndnds,1) = edga(nghO(nbndnds,end),1);
%set first boundary nodes
nghOb_ = nghO(bndnds,:)';
nghOb_(nghOb_==0) = size(edga_,1); edga1_ = edga_(:,1); edga2_ = edga_(:,2);
pot_ = [edga1_(nghOb_); edga2_(nghOb_)]; left = [true(size(nghOb_)); false(size(nghOb_))]; 
C=repmat((1:size(pot_,2)),size(pot_,1),1); [pot_,R] = sort(pot_,1);
I = and(and(and([true(1,size(pot_,2)); pot_(1:end-1,:)~=pot_(2:end,:)],[pot_(1:end-1,:)~=pot_(2:end,:); true(1,size(pot_,2))]),pot_~=0), left(R+(C-1)*size(pot_,1)));
ngh(bndnds,1) = pot_(I);
%count max
nN = sum(nghO~=0,2); nN(bndnds) = nN(bndnds)+1;
R=(1:nds)';
for i=2:maxngh
	C=min(i-1,nN-1);
	ndsrch = ngh(R+(C-1)*size(ngh,1)); ndsrch=reshape(repmat(ndsrch,1,size(nghO,2))',size(nghO_,1),1);
	ngh(R,i) = nghO_(nghO_(:,1) == ndsrch,2).*(nN>=i);
end;