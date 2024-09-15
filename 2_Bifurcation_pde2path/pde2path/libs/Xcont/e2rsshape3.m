function ids=e2rsshape3(p,sig) 
% e2rsshape3: root mean squared (RMS) 
V=p.X; F=p.tri; 
e=edge_lengths(V,F); A=doublearea(V,F)/2; lrms2=sum(e.^2,2)/3;
q=A./lrms2; [qs,idx]=sort(q,'ascend'); fprintf('q_min=%g\n',qs(1)); 
nt=size(p.tri,1); ntr=round(sig*nt); ids=idx(1:ntr); 