function plotbra(varargin)
% PLOTBRA : plot branch of p-structure.
%
%    plotbra(structure,varargin)
% or plotbra(dir,varargin)
% or plotbra(dir,pt,varargin)
%
%    varargin = window number, component number (see below), [option value pairs]
% or varargin = [option value pairs]
% 
% Component number cmp: if number: y-axis, if [x y]: x-, y-axis. 
% counting from the non-user data of bradat(p) onward,e.g. cmp=0 is last 
% component of bradat(p), and cmp entries can be negative.
%
% Option value pairs (option : Description (standard value)): 
%
% *      'lwun'  : Linewidth of unstable sol. (2)
% *      'lwst'  : Linewidth of stable sol. (4)
% *      'lw'    : overrides lwun/lwst: lwst=lw, lwun=lw/2
% *      'tyun'  : Linetype of unstable sol. ('-')
% *      'tyst'  : Linetype of stable Sol. ('-')
% *      'ms'    : Markersize (branch/hopf points) (5)
% *      'rms'   : Markersize (reg. points) (0)
% *      'fms'   : Markersize (fold points) (5)
% *      'lms'   : Markersize (labeled points); if vector [reg., special] markersizes, 
%                  if number markersize for all labeled points ([5,4])
% *      'lsw'   : Global switch to enable/disable label of FP/HP/BP/lab/usrlam; 
%                  e.g. 31 plots all labels
% *      'lab'   : Place markers and labels for given label-list 
%                  ([] if call is plotbra(dir,pt), else all points in dir)
% *      'labi'  : As 'lab', but step through all labels in branch with given increment
% *      'labu'  : boolean for labeling user lambda (false)
% *      'fs'    : Fontsize for labels and axesfs (14)
% *      'lfs'   : Override fontsize for labels, if zero no labels (only markers) (14)
% *      'fp'    : First point to be considered (1)
% *      'lp'    : Last point to be considered (last point in p.branch)
% *      'cl'    : Color, see option in plot ('black')
% *      'bplab' : Place markers and labels for label-list of branch points (all in dir)
% *      'fplab' : Place markers and labels for label-list of fold points (all in dir)
% *      'hplab' : Place markers and labels for label-list of hopf points (all in dir)
% *      'fancy' : 0: no annotation arrow, 1: static annotation arrow, 
%                  2: annotation arrow movable by mouse (1)
% *      'wnr'   : window number (3)
% *      'cmp'   : component number, see above (0)
% *      'os'    : length of the annotation arrow (1)
% *      'odr'   : rough direction of the annotation arrow for regular points ([-1 1])
% *      'ods'   : rough direction of the annotation arrow for special points ([1 0.1])
%
% deprecated:
% *      'inc'   : Show point markers stepping through all points with given increment (1)
%
% See also stanparam, plot, plotbraf (obsolete), bradat.

% set standard parameters
brl=6; % replaces length bradat
wnr=3; cmp=0; xcmp=4; lwun=2; lwst=4; tyun='-'; tyst='-'; 
ms=5; rms=0; fms=5; lms=5; lsms=4; inc=1; lab=[]; labsw=0;
labu=0; % 1: userlam should be labeled
fs=14; lfs=14; fp=1; fp2=0; cl='k'; fancy=1; % annotation style
os=1; % length of x resp y annotation arrows
odr=[-1,1]; % direction of x resp y ann. arrow for regular points. If 0, then random 
ods=[1,0.1]; % like odr for bif/fold/hopf points
nlb=0; nlf=0; nlh=0;  % lengths of labellist of branch points/fold points/hopf points
bplab=[]; fplab=[]; hplab=[]; 
noa=nargin; % number of input-arguments
k=2; % first optional-input
if noa==1;  varargin{2}={'no options'}; end % enable query for varargin{2} if no second input is given.
p=varargin{1}; % check if varargin{1} is struc or dir and load p if necessary
if ischar(p) % call like plotbra('h',...)
   pdir=p; pt=char(['pt' mat2str(max(getlabs(p)))]);
   if ischar(varargin{2})
      if max(strncmp(varargin{2},'pt',2)) || ... % call like plotbra('h','pt5',...)
              max(strncmp(varargin{2},{'fpt','bpt', 'hpt'},3))  
          pt=varargin{2}; k=k+1;  % with pt*, hence no more labels (unless desired)
      end
   end   
   try; p=loadpp(pdir,pt); catch; return; end
end
if ~isstruct(p); error('First input argument must be a (problem) struct or a problem dir.'); end
lsw=1; % default; 
if isfield(p.plot,'lsw'); lsw=p.plot.lsw; end % default, or use entry in p.plot 
% set parameters defined by p
lp=size(p.branch,2); lp2=max(getlabs(p.file.dir));
if isfield(p.plot,'brafig'); wnr=p.plot.brafig; end
if isfield(p.plot,'bpcmp'); cmp=p.plot.bpcmp; end
if isfield(p.plot,'fancybd'); fancy=p.plot.fancybd; end; if fancy==2; fancy=1; end 
if isfield(p.plot,'fs'); fs=p.plot.fs; lfs=p.plot.fs; end
if (k<noa && isnumeric(varargin{k})); % if first inputs are no option pairs, first two inputs are window and component
     wnr=varargin{k}; cmp=varargin{k+1}; k=k+2;
     if length(cmp)>1; xcmp=brl+cmp(1); cmp=cmp(2);  end
end
while k<=noa % parse all optional inputs
    % if input is cell, take the entries as options
    if iscell(varargin{k});
        noa=nargin+length(varargin{k})-1; l=length(varargin)+1; inp=varargin{k};
        for i=k+1:nargin; varargin{l+i-k}=varargin{i}; end
        for i=1:length(inp); varargin{k+i-1}=inp{i}; end
    end
    if ~ischar(varargin{k}); warning(['The input number ' num2str(k) ' could not be interpreted and is ignored. Use "''option'',value" syntax.']);
    end        
    switch lower(varargin{k})
      case 'lwun'; lwun=varargin{k+1}; k=k+2; % Linewidth of unstable Sol.
      case 'lwst'; lwst=varargin{k+1}; k=k+2; % Linewidth of stable Sol.
      case 'lw';   lwst=varargin{k+1}; lwun=lwst/2; k=k+2;
      case 'tyun'; tyun=varargin{k+1}; k=k+2; % Linetype of unstable Sol.
      case 'tyst'; tyst=varargin{k+1}; k=k+2; % Linetype of stable Sol.
      case 'ms';   ms=varargin{k+1}; k=k+2; % Markersize (branch points)           
      case 'rms';  rms=varargin{k+1}; k=k+2; % Markersize (reg. points) 
      case 'fms';  fms=varargin{k+1}; k=k+2; % Markersize (fold points)           
      case 'lms';  lms=varargin{k+1};
                   if length(lms)>1; lsms=lms(2); lms=lms(1);
                   else lsms=lms; end
                   k=k+2; % Markersize (labels) 
      case 'inc';  inc=varargin{k+1}; k=k+2; % increment for point markers
      case 'lab';  lab=varargin{k+1}; labsw=1; k=k+2; % label-list 
      case 'bplab';  bplab=varargin{k+1}; nlb=length(bplab); k=k+2; % branch point label-list, if [] then all
      case 'fplab';  fplab=varargin{k+1}; nlf=length(fplab); k=k+2; % fold point label-list, if [] then all
      case 'hplab';  hplab=varargin{k+1}; nlh=length(hplab); k=k+2; % Hopf point label-list, if [] then all
      case 'labi'; bralen=size(p.branch,2); % labels with increment 
                   lab=1:varargin{k+1}:bralen; lab=p.branch(1,lab); labsw=1; k=k+2; 
      case 'fs';   fs=varargin{k+1}; k=k+2; % fontsize for labels 
      case 'lfs';  lfs=varargin{k+1}; k=k+2; % label fontsize for labels, if 0 no text
      case 'fp';   %fp=varargin{k+1}; k=k+2; % first point
                   fp2=varargin{k+1};fp=min(find(p.branch(1,:)>=varargin{k+1},size(p.branch(1,:),2)));k=k+2;
      case 'lp';   %lp=varargin{k+1}; k=k+2; % last point                   
                   lp2=varargin{k+1}; lp=max(find(p.branch(1,:)<=varargin{k+1},size(p.branch(1,:),2)));k=k+2;
      case 'cl';   cl=varargin{k+1}; k=k+2; % color
      case 'fancy';  fancy=varargin{k+1}; k=k+2; % fancyness
      case 'wnr'; wnr=varargin{k+1}; k=k+2; % window number
      case 'cmp'; cmp=varargin{k+1}; k=k+2; % branch-component(s) to plot
          if length(cmp)>1;  xcmp=brl+cmp(1); cmp=cmp(2); end
      case 'os'; os=varargin{k+1}; k=k+2;
      case 'odr'; odr=varargin{k+1}./max(abs(varargin{k+1}),1e-50); k=k+2;
      case 'ods'; ods=varargin{k+1}./max(abs(varargin{k+1}),1e-50); k=k+2;
      case 'labu'; labu=varargin{k+1}; labusw=1; k=k+2; % labels usrlam if true
      case 'lsw'; lsw=varargin{k+1}; k=k+2;
      case 'auxdict'; p.plot.auxdict=varargin{k+1}; k=k+2; % sets dictionary for components in p.fuha.outfu
      otherwise;
          if ischar(varargin{k});
              warning(['Option ' varargin{k} ' does not exist and is ignored.']);
          end
          k=k+1;
    end
end
if isempty(lab); labu=1; end % lsw is recessive, i.e., 'lab',* dominates lsw=1+x
if fix(lsw/16); lsw=lsw-16; 
    if ~labsw; lab=getlabs(p.file.dir); end % label all regular points 
end 
if fix(lsw/8); lsw=lsw-8; nlf=-1; end % label all folds (if not overwritten by 'fplab',*)
if fix(lsw/4); lsw=lsw-4; nlh=-1; end % label all Hopf 
if fix(lsw/2); lsw=lsw-2; nlb=-1; end % label all BP 
if (lsw>=1 && labu);  lab=[lab,p.branch(1,ismember(p.branch(4,:),p.usrlam))];
   %idx=find(p.branch(2,:)==-3); lab=[lab; p.branch(1,idx)]; 
end % add userlam to label list if desired
lab=unique(lab); nl=length(lab);  % remove multiple entries, save the length of the label-list
pp=p.branch(brl+cmp,:); % set y-Data
figure(wnr); hold on; % set figure and axis
auxdictl=0; if(isfield(p.plot,'auxdict')); auxdictl=length(p.plot.auxdict); end % storage for name of aux para
set(gca,'FontSize',fs);
% xlabel
if(xcmp==4); 
    if(p.nc.ilam(1)<=auxdictl); xax=p.plot.auxdict{p.nc.ilam(1)};
    else xax=['prim.par.=aux.var. ' mat2str(p.nc.ilam(1))]; end
else
    if(xcmp-brl<=auxdictl); xax=p.plot.auxdict{xcmp-brl};  
    else xax=mat2str(xcmp-brl); xax=['user branch comp. ' xax]; end
end
xlabel(char(xax)); 
% ylabel
yax=mat2str(cmp);
if(cmp>0); 
  if(cmp<=auxdictl); yax=p.plot.auxdict{cmp};  
  else yax=['user branch comp. ' yax]; end
else
  switch cmp
    case 0; yax='L2-norm'; case -1; yax='err';
    otherwise; yax=brl+cmp; yax=['branch comp. ' yax];
  end
end
ylabel(char(yax));

% plot branch
for i=fp:lp-1 % loop from first-point to last point
    % linewidth is type of first point, unless this is a fold (2),
    % initial point (-1) (also from swibra (-2)) or branch point (1) or 
    % Hopf (3). In these cases line width is from second point.
    lw=lwun; lt=tyun; 
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=lwst; lt=tyst; end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=lwst; lt=tyst; end
    plot(p.branch(xcmp,i:i+1),pp(i:i+1),lt,'Linewidth',lw,'color',cl); 
end


% mark points
if(rms>0) % plot markers for regular points 
  for i=fp:inc:lp
    if (p.branch(3,i)==0); plot(p.branch(xcmp,i),pp(i),'*','MarkerSize',rms,'color',cl);
    else plot(p.branch(xcmp,i),pp(i),'+','MarkerSize',rms,'color',cl); % unstable point
    end
  end
end

for i=fp:lp; % plot markers for branch/fold/hopf points
    if(ms>0 && ismember(p.branch(2,i),[-2 1])) % plot markers for branch points
      plot(p.branch(xcmp,i),p.branch(brl+cmp,i),'o','MarkerSize',ms*lwun,'color',cl); 
    end
    if(ms>0 && ismember(p.branch(2,i),3)) % plot markers for hopf points
      plot(p.branch(xcmp,i),p.branch(brl+cmp,i),'d','MarkerSize',ms*lwun,'color',cl); 
    end
    if(fms>0 && ismember(p.branch(2,i),2)) % plot markers for fold points
      plot(p.branch(xcmp,i),p.branch(brl+cmp,i),'x','MarkerSize',fms*lwun,'color',cl); 
    end
end

% get axis data
axis tight; h=gca; xlim=get(h,'XLim'); ylim=get(h,'YLim'); 
xmi=xlim(1); xma=xlim(2); ymi=ylim(1); yma=ylim(2); xext=xma-xmi; yext=yma-ymi;

% label regular points
idx=0.1*ones(1,nl); if isempty(lp2); lp2=1e6; end; 
for k=1:nl;
    if (lab(k)>=fp2 && lab(k)<=lp2);
        a=find(p.branch(1,:)==lab(k)); 
        if(~isempty(a)); try; idx(k)=a; catch; end
        elseif p.sw.verb>2
            las=mat2str(lab(k));
            fprintf(['point ' las ' is in directory but not on branch.\n']);
        end
    end
end
idx=idx(idx~=0.1); 
if (~isempty(idx));
    rv=rand(5,length(idx));
    ox=0.02*os*xext*(odr(1)*(1+rv(1,:))+(odr(1)==0)*real(exp(1i*2*pi*rv(2,:))).*(1+rv(3,:))); 
    oy=0.02*os*yext*(odr(2)*(1+rv(4,:))+(odr(2)==0)*imag(exp(1i*2*pi*rv(2,:))).*(1+rv(5,:)));
    x=p.branch(xcmp,idx);y=pp(idx); x1=x+ox;y1=y+oy;  
    xmi=min([x1,xmi]); xma=max([x1,xma]); 
    ymi=min([y1,ymi]); yma=max([y1,yma]);
    xext=xma-xmi;  yext=yma-ymi;
    axis([xmi,xma,ymi,yma]); box on;
    if (fancy==2); sp=get(gca,'Position'); set(gca,'Position',[0 0 1 1]); end 
    for k=1:length(idx); 
        ls=mat2str(p.branch(1,idx(k))); 
        if(lfs>0 && fancy<2); 
        text(x1(k),y1(k),ls,'FontSize',lfs,'color',cl,'VerticalAlignment','bottom'); end;
        if fancy==1; plot([x1(k),x(k)],[y1(k),y(k)],'color',cl); end
        if fancy==2;
              huafixed('textarrow',[x1(k)-0.01*ox(k),x(k),xmi,xext],[y1(k)-0.01*oy(k),y(k),ymi,yext],'string',ls,'FontSize',lfs,'color',cl);
        end
        plot(x(k),y(k),'.','MarkerSize',3*lms*lwun,'color',cl);
    end
    if (fancy==2); set(gca,'Position',sp); end
end
% label branch/fold/hopf points
%bplab, fplab, hplab
if (~(nlb==0) || ~isempty(bplab)); labplot(p,1,bplab,fp2,lp2,cmp,xcmp,cl,lfs,lwun*lsms,fancy,xmi,xma,xext,ymi,yma,yext,os,ods); end % branch pt labels
if (~(nlf==0) || ~isempty(fplab)); labplot(p,2,fplab,fp2,lp2,cmp,xcmp,cl,lfs,lwun*lsms,fancy,xmi,xma,xext,ymi,yma,yext,os,ods); end % fold labels
if (~(nlh==0) || ~isempty(hplab)); labplot(p,3,hplab,fp2,lp2,cmp,xcmp,cl,lfs,lwun*lsms,fancy,xmi,xma,xext,ymi,yma,yext,os,ods); end % hopf labels
end

% function for label of branch/fold/hopf points

function labplot(p,ty,plab,fp,lp,cmp,xcmp,cl,lfs,lwun,...
    fancy,xmi,xma,xext,ymi,yma,yext,os,ods)
% Plot labels for points of type ty from list plab.
% Labels are taken from files in directory to get them right!
% Therefore this function is not used for regular points.
switch ty
    case 1; labs='BP'; lsym='o'; flabs='bpt'; va='top'; 
    case 2; labs='FP'; lsym='x'; flabs='fpt'; va='bottom'; 
    case 3; labs='HP'; lsym='d'; flabs='hpt'; va='bottom';  
end
pts=p.branch(:,p.branch(2,:)==ty);
dn=[p.file.dir '/' flabs '*.mat']; fnames=dir(dn); % lab=[];
fplab=zeros(1,length(fnames));
for j=1:length(fnames) % create list of branch point labels
    cn=fnames(j).name; fplab(j)=str2double(cn(4:length(cn)-4)); 
end
if (isempty(plab)); plab=fplab; end % if empty then take all
% test if plab entries are valid
x=0.5*ones(1,length(plab));
existing=ones(1,length(plab));
for k=1:length(x);
    try x(k)=pts(xcmp,plab(k));
    catch; existing(k)=0; if p.sw.verb>=2; fprintf([flabs num2str(plab(k)) ' does not exist.\n']); end
    end
end
x=x(x~=0.5);
existing=logical(existing);
plab=plab(existing);
nl=length(plab);
if(isempty(plab)); return; end % if still empty, e.g. there are no branch/fold/hopf points stop
% plot branch/fold/hopf labels
rv=rand(5,nl); % 5 random vectors rand(1,nl)
ox=0.05*os*xext*(ods(1)*(1+rv(1,:))+(ods(1)==0)*real(exp(1i*2*pi*rv(2,:))).*(1+rv(3,:))); oy=0.1*os*yext*(ods(2)*(1+rv(4,:))+(ods(2)==0)*imag(exp(1i*2*pi*rv(2,:))).*(1+rv(5,:)));
pp=pts(brl+cmp,:); % the user data
y=pp(plab);

x1=x+ox; y1=y+oy; xmi=min([x1,xmi]); xma=max([x1,xma]);
ymi=min([y1,ymi]); yma=max([y1,yma]); xext=xma-xmi; yext=yma-ymi;
axis([xmi,xma,ymi,yma]);
if (fancy==2); sp=get(gca,'Position'); set(gca,'Position',[0 0 1 1]); end 
for k=1:nl; 
    if (pts(1,k)>=fp && pts(1,k)<=lp && ismember(plab(k),fplab)) 
      ls=[labs mat2str(plab(k))];
      if(lfs>0) 
         switch fancy
             case 0;  text(x(k),y(k),ls,'FontSize',lfs,'color',cl,...
                'horizontalalignment','center','VerticalAlignment',va); 
             case 1; % with lines, not movable by mouse
                text(x1(k),y1(k),ls,'FontSize',lfs,'color',cl,...
                'horizontalalignment','center','VerticalAlignment',va); 
                plot([x1(k),x(k)],[y1(k),y(k)],'color',cl); 
             case 2; % via annotation-fix, somewhat movable by mouse
                huafixed('textarrow',[x1(k)-0.01*ox(k),x(k),xmi,xext],[y1(k)-0.01*oy(k),y(k),ymi,yext],'string',ls,'FontSize',lfs,'color',cl);
         end         
      end;
      plot(x(k),y(k),lsym,'MarkerSize',lwun,'color',cl);
    end
end
if (fancy==2); set(gca,'Position',sp); end
end