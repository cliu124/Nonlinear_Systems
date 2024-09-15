function hotanplot(p,wnr,cmp,pstyle)
% hotanplot: plot Hopf tangent
% hotanplot(p,wnr,cmp,pstyle)
p.hopf.y=reshape(p.hopf.tau(1:end-2),p.nu,p.hopf.tl); hoplot(p,wnr,cmp,pstyle); 
p.hopf.y=p.mat.M\p.hopf.y0d; hoplot(p,wnr+1,cmp,pstyle); 