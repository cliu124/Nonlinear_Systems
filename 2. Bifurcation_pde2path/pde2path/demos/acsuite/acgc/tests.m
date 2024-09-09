%% check if SM and full jac yield the same EVals
p=loadp('1D1','pt10'); %p=loadp('2Dh8','pt30'); 
%% full 
p.sw.eigssol=0; p.jfac=1; r=resi(p,p.u); Gu=getGu(p,p.u,r); p.sw.verb=2; 
tic; [ineg,muv1,V]=spcalc(Gu,p); toc; muv1
%% SM
p.sw.eigssol=1; p.jfac=0; r=resi(p,p.u); Gu=getGu(p,p.u,r); p.sw.verb=2; 
tic; [ineg,muv2,V]=spcalc(Gu,p); toc; muv2