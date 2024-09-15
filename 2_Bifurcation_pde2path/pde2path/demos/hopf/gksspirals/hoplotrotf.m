function hoplotrotf(dir,sfname,wnr,cnr,pstyle)
try p=loadp(dir,sfname); catch; fprintf(['hoplotrotf:' dir '/' sfname ' does not exist.\n']); return; end  
fprintf('lam=%g\n',getlam(p)); 
%hoplotrot3(p,wnr,cnr,pstyle);  
hoplotrot4(p,wnr,cnr,pstyle);  