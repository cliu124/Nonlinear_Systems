function pmove(odir,ndir)
% PMOVE: move p2p dir data and change file name variables
%
%  pmove(odir,ndir)
%
% See also cp2pdir, setfn
ok=pcopy(odir,ndir);
if(ok==1) 
    rmdir(odir, 's');
end
