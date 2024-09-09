function p=resetc(p)
% RESETC: reset counters in p and clear p.branch
%
%  p=resetc(p)
p.branch=[]; p.file.count=0; p.file.bcount=1; p.file.fcount=1; p.file.hcount=1;
