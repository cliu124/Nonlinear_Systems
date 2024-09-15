function vl=vallist(p) 
% vallist: list of node valences
A=adjacency_matrix(p.tri); 
vl=sum(A,2); vl=(full(vl))'; 