function p=retrigX(p) 
% retrigX: re-triangulate X, taken from distmesh 
t=p.tri; [t2t,t2n]=mkt2t(t); %  element connectivities
t2t=int32(t2t-1)'; t2n=int8(t2n-1)';
[t,t2t,t2n]=trisurfupd(int32(t-1)',t2t,t2n,p.X');  % Update triangles (mex) 
t=double(t+1)'; p.tri=t; p.nt=size(t,1); 
