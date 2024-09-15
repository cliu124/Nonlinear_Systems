function plotsolf(dir,sfname,wnr,cnr,pstyle)
%PLOTSOLF: plot component of p.u from file dir/sfname.mat 
%  plotsolf(dir,sfname,wnr,cnr,pstyle)
% cnr=component number, wnr=window number, pstyle=plot style
%
% See also plotsolu, plotsol
p=loadp(dir,sfname); fprintf('lam=%g\n',getlam(p)); 
plotsol(p,wnr,cnr,pstyle);  tname=[dir '/' sfname]; 
title(['P at ' tname]); end 
