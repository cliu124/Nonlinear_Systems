function p=loadp3(p,dir,fname)
% LOADP3: overwrite p.u from different file
%
%  p=loadp3(p,dir,fname)
%
% loads file dir/fname.mat  which only contains u (and tint time data t,ts) 
% See also loadp.
ffname=[dir '/' fname '.mat'];
s=load(ffname,'p'); q=s.p; p.t=q.t; p.ts=q.ts; p.u=q.u; 
