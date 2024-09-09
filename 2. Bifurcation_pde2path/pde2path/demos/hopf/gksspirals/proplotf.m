% "profile plot" from file, i.e., plot of first time-slice
function proplotf(dir,sfname,wnr,cnr,pstyle)
p=loadp(dir,sfname); fprintf('lam=%g, T=%g\n',getlam(p), p.hopf.T); 
proplot(p,wnr,cnr,pstyle);  
tname=[dir '/' sfname ', T=' mat2str(p.hopf.T,3)]; title(tname);
