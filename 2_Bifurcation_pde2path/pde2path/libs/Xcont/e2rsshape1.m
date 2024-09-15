function ids=e2rsshape1(p,sig) 
% e2rsshape4: list of r/h (incircle/max-edge-length) 
V=p.X; F=p.tri; e=edge_lengths(V,F); el=max(e,[],2); % long edges 
s=sum(e,2)/2; A=doublearea(V,F)/2; r=A./s; 
qi=r./el; [qis,idx]=sort(qi,'ascend');
nt=size(p.tri,1); ntr=round(sig*nt); ids=idx(1:ntr); 

