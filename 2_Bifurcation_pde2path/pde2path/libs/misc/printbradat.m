function printbradat(pdir,varargin)
% PRINTBRADAT: print data from branch of file or p
%
%  printbradat(p)        - all of p.branch
%  printbradat(dir)      - use dir/ptN with largest N
%  printbradat(dir,pt)   - use dir/pt
%  printbradat(dir,pt,ilab) - ilab>0 only label of p in dir/pt
%  (ilab=2: do not plot header, ilab=3 only plot data)
%
% See also printaux, getaux, bradat, stanbra
ilab=0; lab=-1;
if ischar(pdir); 
    if isempty(varargin); pt=char(['pt' mat2str(max(getlabs(pdir)))]); 
    else pt=varargin{1}; end
    try p=loadp(pdir,pt); 
    catch fprintf(['point ' pdir '/' pt ' does not exist.\n']); 
        pt=ptselect(pdir); p=loadp(pdir,pt); 
    end
    if nargin>=3; 
        ilab=varargin{2};
        lab=p.file.count-1;
        if(ilab==1) 
            fprintf(['Branch data of ' pdir '/' pt]); 
            fprintf([' (label ' mat2str(lab) ' only)\n']);
        end
    else
        fprintf(['Branch data of ' pdir '/' pt '\n']); 
    end
else
    p=pdir; fprintf(['Branch data of ' inputname(1) '\n']);
end
auxv=[]; auxf=[]; auxdictl=0;
if (isfield(p.plot,'auxdict')) auxdictl=length(p.plot.auxdict); end
for i=1:length(p.u(p.nu+1:end));
    if(i<=auxdictl)
        auxv=[auxv p.plot.auxdict{i} blanks(11-length(p.plot.auxdict{i}))];
    else
        auxv=[auxv 'aux' mat2str(i) '       '];
    end
    auxf=[auxf '%5.2e '];
end
if(auxdictl>=p.nc.ilam(1)) 
    primaux=p.plot.auxdict{p.nc.ilam(1)};
else
    primaux='lambda';
end
if(ilab<3) 
  %  fprintf(['step type ineg  ' primaux blanks(11-length(primaux)) 'err        L2-norm    ' auxv 'max|u_1|   min|u_1|\n']);
  fprintf(['step type ineg \n']); 
end
s=size(p.branch); l=s(2);
for i=1:l;
  if(lab<0 || (lab>=0 && p.branch(1,i)==lab))
    fprintf('%4i %4i %4i',p.branch(1,i),p.branch(2,i),p.branch(3,i));
    for j=4:length(p.branch(:,i));
        fprintf(' % 5.2e ',p.branch(j,i));
    end
    fprintf('\n');
  end
end
if(ilab<2) 
    out=['Point types: -1=initial,-3=usrlam,-2=swibra,0=normal,'...
        ,'1=branch point,2=fold,3=Hopf,4=normal on Hopf branch\n'];
    fprintf(out);
end