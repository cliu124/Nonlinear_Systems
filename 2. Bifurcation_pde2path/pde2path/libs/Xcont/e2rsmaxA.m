function ids=e2rsmaxA(p,sig)
% e2rsmaxA:  ElementToRefineSelector of tri with |tri|>(1-sig)*p.maxA
A=doublearea(p.X,p.tri)/2; ids=find(A>(1-sig)*p.maxA); 