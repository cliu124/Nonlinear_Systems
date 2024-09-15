function p=fchpostmeshmod(p)
%global eta; 
% post mesh-modification procedure
try se=size(p.eta,2); catch; p.eta=[]; se=0; end
if(se~=size(p.u,1)) % eta not yet set, or mesh is refined 
  C=n2triamat(p.mesh.p,p.mesh.t); ta=triar(p.mesh.p,p.mesh.t); 
  eta=ta*C; p.eta=[eta zeros(1,p.np)]; 
end
if(p.sw.sfem~=0) p=setfemops(p); end 
end 
