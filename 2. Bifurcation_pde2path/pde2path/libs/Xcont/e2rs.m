function ids=e2rs(p,sig)
% e2rs: ElementToRefineSelector based on area 
A=doublearea(p.X,p.tri)/2; [As,idx]=sort(A,'descend'); 
nt=size(p.tri,1); ntr=round(sig*nt); ids=idx(1:ntr); 