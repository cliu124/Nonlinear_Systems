function [ineg,muv]=vspcalc(Gu,p)
% vspcalc: loop over shifts and call spcalc 
%
%  [ineg,muv]=vspcalc(Gu,p)
for j=1:length(p.nc.eigref); [ineg(j),muv(j,:),V]=spcalc(Gu,p,j); end
