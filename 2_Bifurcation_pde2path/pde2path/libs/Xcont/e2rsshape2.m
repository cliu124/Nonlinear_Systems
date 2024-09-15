function ids=e2rsshape2(p,sig) 
% e2rsshape2: list of smallest A/(Sum(edgelengths))^2/3 
V=p.X; F=p.tri; 
e=edge_lengths(V,F);A=doublearea(V,F); s=sum(e,2);
q=A./(s.^(2/3)); [qs,idx]=sort(q,'ascend'); fprintf('q_min=%g\n',qs(1)); 
nt=size(p.tri,1); ntr=round(sig*nt); ids=idx(1:ntr); 

