function p=doubledrop(p)
% doubledrop: double the size of p.mat.drop and p.mat.fill, e.g. for BP continuation 
r=size(p.mat.drop,1); s=size(p.mat.drop,2);
if r>1;  p.mat.drop = [[p.mat.drop zeros(r,s)]; [zeros(r,s) p.mat.drop]]; 
  p.mat.fill = [[p.mat.fill zeros(s,r)]; [zeros(s,r) p.mat.fill]]; 
end