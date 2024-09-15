function plot1Df(dir,pt,wnr,wfak,efak,hfak,sw)
p=loadp(dir,pt); fprintf('lam=%g\n',getlam(p)); 
tit=[dir '/' pt]; plot1D(p,wnr,wfak,efak,hfak,sw,tit); 
valf(dir,pt)
