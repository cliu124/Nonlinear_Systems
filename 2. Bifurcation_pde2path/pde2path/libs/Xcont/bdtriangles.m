function idx=bdtriangles(T)
% bdtriangles: find boundary triangle indices 
F=boundary_faces(T); 
[l1,idx]=ismember(F,T(:,[1 2]),'rows'); idx1=idx(l1); 
[l2,idx]=ismember(F,T(:,[2 3]),'rows'); idx2=idx(l2); 
[l3,idx]=ismember(F,T(:,[3 1]),'rows'); idx3=idx(l3); 
idx=unique([idx1; idx2; idx3;]); 