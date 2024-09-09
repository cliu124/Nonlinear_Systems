function p=loadp2(dir,fname,fname0,varargin)
% LOADP2: load p.u, and other data (e.g., mesh) from different files,
%
%  p=loadp2(dir,fname,fname0,varargin)
%
% loads file dir/fname.mat  which only contains u (and tint time data t,ts) 
% other data from dir/fname0.mat
%
% See also loadp.
noa=nargin-3; % number of opt arguments 
p=loadp(dir,fname0);
ffname=[dir '/' fname '.mat'];
s=load(ffname,'p'); q=s.p; p.t=q.t; p.ts=q.ts; p.u=q.u; 
if noa>0; p=setfn(p,varargin{1}); end;