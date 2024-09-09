function [VV,FF]=rmdegfaces(V,F,varargin)
% rmdegfaces: mod of REMOVE_DEGENERATE_FACES from gptoolbox 
  max_iter=100; epsilon=0; keepbd=0; 
  params_to_variables=containers.Map( ...
    {'MaxIter','Epsilon','keepbd'}, {'max_iter','epsilon','keepbd'}); %HU 
  v=1;
  while v <= numel(varargin)
    param_name=varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin)); v=v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else; error('Unsupported parameter: %s',varargin{v}); end
    v=v+1;
  end  
  [bdF,bdJ,bdK]=boundary_faces(F); [V,FF,~,J]=faces_first(V,F,bdF); EJ=(1:size(V,1))';
  % function for Combinatorially degenerate faces
  combinatorially_degenerate=@(G) any(G(:,2)==G(:,3) | G(:,3)==G(:,1) | G(:,1)==G(:,2),2);
 % if keepbd; D(bdtrili(FF))=0; end 
  iter=1;
  get_edge=@(G,I,S) ...
      G([sub2ind(size(G),find(I),mod(S+1-1,3)+1) sub2ind(size(G),find(I),mod(S+2-1,3)+1)]);
  validate_edges= @(F,E) assert( ...
    all(ismember(sort(E,2),sort([F(:,[2 3]);F(:,[3 1]);F(:,[1 2])],2),'rows')));
  while true
    [IA,~,l,A]=is_acute(V,FF); 
    [~,bdJ,~]=boundary_faces(FF); 
    D=A/2 < epsilon; 
    if keepbd; D(bdtrili(FF))=0; end 
    % Find non-delaunay, longest edges of small, non-acute triangles
    L=~is_intrinsic_delaunay(V,FF) & bsxfun(@eq,l,max(l,[],2)) & repmat(D&~IA,1,3); 
    % bsxfun: apply 'equal' element-wise 
    [LI,LJ]=ind2sub(size(FF),find(L));
    EL=FF([sub2ind(size(FF),LI,mod(LJ+1-1,3)+1) sub2ind(size(FF),LI,mod(LJ+2-1,3)+1) ]);
    % Include all "non flippable" obtuse cases in collapsable
    IC=D & ~any(L,2);  % is collapseble? 
    %if keepbd; IC(bdtrili(FF))=0; end 
    % Find short edges of acute triangles
    [~,S]=min(l(D&IC,:),[],2); ES=get_edge(FF,D&IC,S);

    if size(ES,1) == 0 && size(EL,1) == 0; break; end
    validate_edges(FF,ES); validate_edges(FF,EL);
    [ESL,EI]=conservative_edge_matching([ES;EL],'Method','recursive');    % 
    EL=ESL(EI>size(ES,1),:); ES=ESL(EI<=size(ES,1),:);
    validate_edges(FF,ES);   validate_edges(FF,EL);
    % Flipping an edge can never create a degenerate face so do those first
    validate_edges(FF,EL);   FF=flip_edges(FF,EL);
    validate_edges(FF,ES);   ES=sort(ES,2);
    EJ(ES(:,2))=ES(:,1);
    %if keepbd; EJ=rmbdtrifromlist(tri,ids) EJ=union(EJ,
    FF=EJ(FF);   % the actual reduction! 
    FF=FF(~combinatorially_degenerate(FF),:); %size(FF), pause 

    if iter == max_iter; warning('Maximum iterations reached'); break; end
    iter=iter + 1;
  end

  [VV,UJ,I]=remove_unreferenced(V,FF); 
  FF=UJ(FF); J=UJ(EJ(J));

end
