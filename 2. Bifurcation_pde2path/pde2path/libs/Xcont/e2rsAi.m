function ids=e2rsAi(p,sig)
% e2rsAi:  selecting small areas (mainly for coarsening) 
A=doublearea(p.X,p.tri)/2; [~,idx]=sort(A,'ascend'); 
nt=size(p.tri,1); ntr=round(sig*nt); ids=idx(1:ntr); 