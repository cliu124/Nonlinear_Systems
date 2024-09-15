function printaux(p,varargin)
% PRINTAUX: print auxiliary variable index and value
%
%  printaux(p)      - all aux vars of p
%  printaux(dir)    - all aux vars of dir/ptN with largest N
%  printaux(dir,pt) - all aux vars of dir/pt
%  printaux(p,1)    - only active aux vars of p
%
% if varargin nonempty print only primary and other active ones
act=0; % act=active=0 for all, 1 for only active ones 
if ischar(p) 
    pdir=p;
    if isempty(varargin); pt=char(['pt' mat2str(max(getlabs(pdir)))]); 
    else pt=varargin{1}; end
    try p=loadp(pdir,pt); 
    catch fprintf(['point ' pdir '/' pt ' does not exist.\n']); 
        pt=ptselect(pdir); p=loadp(pdir,pt); 
    end
    fprintf(['aux vars of ' pdir '/' pt ' (is point of type ' mat2str(p.branch(2,end)) ')\n']);
elseif(isempty(p.branch))
    fprintf('Branch empty.\n'); 
else
    fprintf(['aux vars of ' inputname(1)  ' (is point of type ' mat2str(p.branch(2,end)) ')\n']);
    if ~isempty(varargin); act=1; end   
end
auxdictl=0; if(isfield(p.plot,'auxdict')); auxdictl=length(p.plot.auxdict); end % storage for name of aux para
for k=1:auxdictl, lab{k}=p.plot.auxdict{k}; end
for k=auxdictl+1:length(p.u)-p.nu, lab{k}=mat2str(k); end
if(act==0)
    fprintf('index label  value\n');
    for k=p.nu+1:length(p.u) fprintf(' %i %s   %g \n',k-p.nu,lab{k-p.nu}, p.u(k)); end
else
    fprintf('p.nc.ilam  values\n');
try;  
    for k=1:p.nc.nq+1; 
   fprintf('  %i  %s   %g \n',p.nc.ilam(k),lab{p.nc.ilam(k)},p.u(p.nu+p.nc.ilam(k))); 
  %    fprintf('  %i  %g \n',p.nc.ilam(k),p.u(p.nu+p.nc.ilam(k))); 
    end 
catch; end 
end
end
