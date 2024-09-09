function c3=c3fu(c,nt)
% c3fu: generate 3-node triangulation from 6-node 
% Adapted from Pozrikidis' FSElib, obsolete, see tri2six 
Ic=0; 
for i=1:nt
 Ic=Ic+1; c3(Ic,1)=c(i,1); c3(Ic,2)=c(i,4); c3(Ic,3)=c(i,6);
 Ic=Ic+1; c3(Ic,1)=c(i,4); c3(Ic,2)=c(i,2); c3(Ic,3)=c(i,5);
 Ic=Ic+1; c3(Ic,1)=c(i,5); c3(Ic,2)=c(i,3); c3(Ic,3)=c(i,6);
 Ic=Ic+1; c3(Ic,1)=c(i,4); c3(Ic,2)=c(i,5); c3(Ic,3)=c(i,6);
end