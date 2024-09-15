function ok=pcopy(odir,ndir)
% PCOPY: copy p2p dir data and change file name variables
%
%  pcopy(odir,ndir)
%
% See also mvp2pdir, setfn
ok=1;
% create new dir
if ~exist(odir,'dir') 
    fprintf('Directory %s does not exists.\n',odir);
    ok=0; return;
end
if exist(ndir,'dir')
    fprintf('Directory %s already exists.\n',ndir);
    ok=0; return;
end
fprintf('Warning: copying first level of directory only!\n');
mkdir(ndir); fprintf('creating directory %s\n',ndir);
% loop through files and change file name variables
dn=getlabs(odir); dn=sort(dn);
for j=1:length(dn)
    cn=['pt' mat2str(dn(j))]; 
    chg(odir,cn,ndir);
end
dn=[odir '/bpt*.mat']; c=dir(dn);
for j=1:length(c)
    cn=c(j).name; lab=cn(4:length(cn)-4); cn=['bpt' lab]; 
    chg(odir,cn,ndir);
end
dn=[odir '/fpt*.mat']; c=dir(dn);
for j=1:length(c)
    cn=c(j).name; lab=cn(4:length(cn)-4); cn=['fpt' lab]; 
    chg(odir,cn,ndir);
end
dn=[odir '/hpt*.mat']; c=dir(dn);
for j=1:length(c)
    cn=c(j).name; lab=cn(4:length(cn)-4); cn=['hpt' lab]; 
    chg(odir,cn,ndir);
end
end
% Change file name variables and save
function chg(odir,cn,ndir)
    p=loadpp(odir,cn);
    p.file.dir=ndir; p.file.pname=[ndir '/pt']; 
    p.file.bpname=[ndir '/bpt']; p.file.hpname=[ndir '/hpt'];  
    p.file.fpname=[ndir '/fpt']; 
    p.file.count=p.file.count-1; % correct increase from loadp 
    if(p.sol.ptype==1) p.file.bcount=p.file.bcount-1; end
    if(p.sol.ptype==2) p.file.fcount=p.file.fcount-1; end
    if(p.sol.ptype==3) p.file.hcount=p.file.hcount-1; end
    fname=[ndir '/' cn]; save(fname,'p'); 
end
