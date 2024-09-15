function printdirdat(pdir)
% PRINTDIRDAT: print data for each point stored in directory
% including branch, fold and hopf points
%
%  printdirdat(dir)
%
% See also printbradat

fprintf(['Data of ptN in directory ' pdir '\n']);
dn=getlabs(pdir); dn=sort(dn);
for j=1:length(dn)
    cn=['pt' mat2str(dn(j))]; 
    if(j==1) printbradat(pdir,cn,2); end
    if(j>1) printbradat(pdir,cn,3); end
end
dn=[pdir '/bpt*.mat']; c=dir(dn);
for j=1:length(c)
    cn=c(j).name; lab=cn(4:length(cn)-4); 
    cn=['bpt' lab]; printbradat(pdir,cn,1);
end
dn=[pdir '/fpt*.mat']; c=dir(dn);
for j=1:length(c)
    cn=c(j).name; lab=cn(4:length(cn)-4); 
    cn=['fpt' lab]; printbradat(pdir,cn,1);
end
dn=[pdir '/hpt*.mat']; c=dir(dn);
for j=1:length(c)
    cn=c(j).name; lab=cn(4:length(cn)-4); 
    cn=['hpt' lab]; printbradat(pdir,cn,1);
end
end
