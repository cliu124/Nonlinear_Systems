function plotsolf(dir,sfname,wnr,cnr,pstyle)
p=loadp(dir,sfname); fprintf('lam=%g\n',getlam(p)); 
plotsol(p,wnr,cnr,pstyle); tname=[dir '/' sfname]; xlabel('x'); 
if 0 
%title(['P at ' tname]); 
[h,hp]=hfu(p,p.u); k=cficon(p,p.u); jc=cfijcf(p,p.u); %l=p.u(p.np+1)
tit={[tname '_., k=(' mat2str(k(1),4) ', ' mat2str(k(2),4) ')']; 
    ['h=(' mat2str(h(1),4) ', ' mat2str(h(2),4) '), J_c=' mat2str(jc,4)]}; 
set(gca,'FontSize',p.plot.fs);
title(tit); 
else 
    gtext(tname,'FontSize',14);
end