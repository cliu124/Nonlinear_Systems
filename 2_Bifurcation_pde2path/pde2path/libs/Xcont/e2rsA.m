function ids=e2rsA(p,sig)
% e2rsA:  ElementToRefineSelector based on area 
if size(p.tri,2)==2; M=massmatrix(p.X,p.tri,'voronoi'); A=sum(M,1); 
else A=doublearea(p.X,p.tri)/2; 
end
[~,idx]=sort(A,'descend'); nt=size(p.tri,1); ntr=round(sig*nt); ids=idx(1:ntr); 