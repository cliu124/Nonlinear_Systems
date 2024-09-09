function labs=getlabs(pdir)
% GETLABS: get point label list from dir
%
%  labs=getlabs(pdir) 
%
% See also setfn, plotbraf
dn=[pdir '/pt*.mat'];c=dir(dn);labs=[];
for j=1:length(c)
    cn=c(j).name; lab=cn(3:length(cn)-4); lab=str2num(lab);
    labs=[labs; lab];
end
