%%
%p=loadp('1D1','pt50'); 
p=loadp('2Dh8','pt30'); 
%%
r=resi(p,p.u); Gu=getGu(p,p.u,r); p.sw.verb=2; 
%%
p.sw.eigssol=0; tic; [ineg,muv1,V]=spcalc(Gu,p); toc; muv1
p.sw.eigssol=1; p.jfac=0; tic; [ineg,muv2,V]=spcalc(Gu,p); toc; muv2