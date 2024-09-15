function plot1Df(dir,pt,wnr,ufac,jfac,sw)
p=loadp(dir,pt); fprintf('lam=%g\n',getlam(p)); 
tit=[dir '/' pt]; %tit=[]; 
plot1D(p,wnr,ufac,jfac,sw,tit); 


