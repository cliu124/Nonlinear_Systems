function [z,p]=blssbel(A,w,p)
% BLSSBEL: Increases p.bel.bw and calls lssbel.
% When the border width p.bel.bw>1, it is recommended to set 
% p.fuha.lss=@lssbel and p.fuha.blss=@blssbel.
%
%  [z,p]=lssbel(A,w,p)
%
% see also: lss, blss, lssbelpi, lssbss, bel, bss
p.bel.bw=p.bel.bw+1;
[z,p]=lssbel(A,w,p);
p.bel.bw=p.bel.bw-1;