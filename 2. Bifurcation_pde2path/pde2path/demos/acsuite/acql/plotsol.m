function plotsol(varargin) % pimp up of old plotsol which also replaces 
% plotsol: plot component of p.u in struct p.
%  plotsol(p,varargin) 
% plotsolf if first argument is a string, i.e., calls like 
% plotsol('h','pt4',options) or plotsol('h',options) supported
% see 'plot' tutorial 
if ischar(varargin{1});
    if nargin>1 && ischar(varargin{2}) %check if varargin{2} is a specific point in folder
    dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
    if exist(str,'file')==2
        p=loadp(dir,varargin{2}); anf=3;  %folder and point is given
    else  % check if user tried to load point which does not exist
         % (something like plotsol('h','pt4')), or if user gives only folder and wants to
         % load the last point in the folder (call like plotsol('h','levn',3))
        if (pt(1)=='p' && pt(2)=='t') || (pt(1)=='f' && pt(2)=='p' && pt(3)=='t') || (pt(1)=='b' && pt(2)=='p' && pt(3)=='t')
            try p=loadp(dir,pt); anf=3; 
            catch; fprintf(['point ' dir '/' pt ' does not exist.\n']); 
            npt=ptselect(dir); p=loadp(dir,npt); pt=npt; anf=3; end
        else; sfname=char(['pt' mat2str(max(getlabs(dir)))]);
            p=loadp(dir,sfname);anf=2;
        end
    end
    else %user called something like plotsol('h',3,1,2) or plotsol('h')
            dir=varargin{1}; sfname=char(['pt' mat2str(max(getlabs(dir)))]);
            p=loadp(dir,sfname); anf=2; 
    end % if nargin>1 && ischar(varargin{2})
    fprintf('lam=%g\n',getlam(p)); % show lambda in the work space
else p=varargin{1}; anf=2;% if first entry is a structure
end
%check existence of fields in p.plot and set fig. nr., component nr, and pstyle 
if isfield(p.plot,'pfig'); wnr=p.plot.pfig; else wnr=1; end
if isfield(p.plot,'pcmp'); cnr=p.plot.pcmp; else cnr=1; end
if cnr==0; cnr=1; end % I do not why but sometime p.plot.pcmp=0 for saved sol.
if isfield(p.plot,'pstyle'); pstyle=p.plot.pstyle; else pstyle=1; end
% check if fig. nr., component nr, and pstyle are passed in form plotsol(p,11,1,2)
if nargin>=anf && isa(varargin{anf},'double'); wnr=varargin{anf}; end;
if nargin>=anf+1 && isa(varargin{anf},'double') && isa(varargin{anf+1},'double'); 
    cnr=varargin{anf+1};
end
if nargin>=anf+2 && isa(varargin{anf},'double') && isa(varargin{anf+1},'double')...
        && (isa(varargin{anf+2},'double')) 
    pstyle=varargin{anf+2};
end
% read out general options
for k=1:nargin
 if ischar(varargin{k})
   switch lower(varargin{k})
     case 'pfig'; wnr=varargin{k+1}; case 'pcmp'; cnr=varargin{k+1};
     case 'pstyle'; pstyle=varargin{k+1}; case 'fs'; p.plot.fs=varargin{k+1};
     case 'cm'; p.plot.cm=varargin{k+1}; case 'axis'; p.plot.axis=varargin{k+1};
   end
 end
end
% set fig. nr. and u
figure(wnr); upde=p.mat.fill*p.u(1:p.nu); 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; u=upde(n0:n1); 
clf(wnr);
try; tname=[p.file.pname,num2str(p.file.count-1)]; catch tname=[]; end 
if p.sol.ptype==1; tname=[p.file.bpname,num2str(p.file.bcount-1)]; end
if p.sol.ptype==2; tname=[p.file.fpname,num2str(p.file.fcount-1)]; end
[po,~,~]=getpte(p);
switch size(po,1)
  case 1
    for k=1:nargin  % read out extra 1D-options    
      if ischar(varargin{k})
        switch lower(varargin{k})
          case 'lw'; p.plot.lw=varargin{k+1}; case 'cl'; p.plot.cl=varargin{k+1};
         end
      end
    end
    if isfield(p.plot,'lw')==0; p.plot.lw=1; end
    if isfield(p.plot,'cl')==0; p.plot.cl=zeros(length(cnr),3); end
    psc=cell(1,length(pstyle)); % plot style cell
    if isa(pstyle,'double')
      for i=1:length(pstyle)
        switch pstyle(i)
          case 1; psc{i}='-'; case 2; psc{i}='--'; case 3; psc{i}=':'; 
              case 4; psc{i}='*';
        end
      end
    else % convert pstyle to the cell {pstyle} 
      if isa(pstyle,'cell')==0; psc={pstyle}; else psc=pstyle; end
    end               
    if iscell(p.plot.cl); cl=zeros(length(p.plot.cl),3); % set colors
     for i=1:length(p.plot.cl); cl(i,:)=pde2pathcolors(p.plot.cl{i}); end
        p.plot.cl=cl;
    end
    if ischar(p.plot.cl); p.plot.cl=pde2pathcolors(p.plot.cl); end      
  %%%%%%%%%%%%%%%% begin filling %%%%%%%%%%%%%%%%%%%%%%%
  % if plotsol(p,11,[1,2],1) instead of plotsol(p,11,[1,2],[1 1]), we fill psc
    if length(cnr)>length(psc); dum=cell(1,length(cnr));
       for i=1:length(cnr); dum{i}=psc{1}; end; psc=dum;
    end 
  % if user types plotsol(p,11,[1,2],[1,2],'lw',3) 
  % instead of plotsol(p,11,[1,2],[1 2],'lw',[3,3]), we fill p.plot.lw
    if length(cnr)>length(p.plot.lw); dum=zeros(1,length(cnr));
      for i=1:length(cnr); dum(i)=p.plot.lw(1); end; p.plot.lw=dum;
    end 
  % if user types plotsol(p,11,[1,2],'cl','b3') instead of 
  % plotsol(p,11,[1,2],'lw',{'b3','b3'}), we fill p.plot.lw
    if length(cnr)>size(p.plot.cl,1); cl=zeros(length(cnr),3);
      for i=1:length(cnr); cl(i,:)=p.plot.cl(1,:); end; p.plot.cl=cl;
    end %%%%%%%% end filling, now plot solution(s)
    for i=1:length(cnr);
       n0=(cnr(i)-1)*p.np+1; n1=cnr(i)*p.np; 
       p.pdeo.grid.plot(upde(n0:n1),psc{i},'linewidth',p.plot.lw(i),...
               'color',p.plot.cl(i,:));
    end
    drawnow; str=[]; % create string like u1, u2 
    if p.nc.neq==1; str='u'; % if we have one component system 
    else
      for i=1:length(cnr)-1; 
          str=[str,'\color[rgb]{',num2str(p.plot.cl(i,:)),'}u',num2str(cnr(i)),', ']; 
      end;
      if cnr==1; str=[str,'u',num2str(cnr(end))];  %use color black for u for 1 comp. plots
      else; str=[str,'\color[rgb]{',num2str(p.plot.cl(end,:)),'}u',num2str(cnr(end))]; % use different colors for u1 and u2         
      end
    end
    set(gca,'FontSize',p.plot.fs); title([str,'\color{black} at ',tname]);
    xlabel(''); axis tight; box on; % set title, box, and label     
  case 2; [po,t,e]=getpte(p); %%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%
  % check if oopde or pdetoolbox is used and use resp plotting functions
      if p.sw.sfem==0 || p.sw.sfem==1
        if pstyle==0; pdemesh(po,e,t); end 
        if pstyle==1; pdemesh(po,e,t,u); end 
        if pstyle==2; pdeplot(po,e,t,'xydata',u); end 
        if pstyle==3; h=pdesurf(po,t,u); view(10,40);
        light('Position',p.plot.lpos,'Style','local','Color',[1 1 0]);
        lighting phong; 
        set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8); 
        end        
%   if p.plot.labelsw; xlabel('x','FontSize',p.plot.fs); ylabel('y','FontSize',p.plot.fs); end;
      else
        if pstyle==0; p.pdeo.grid.plot; end
        if pstyle==1; p.pdeo.grid.plot(u,'FaceColor','white'); end 
        if pstyle==2; p.pdeo.grid.plot(u,'LineStyle','none'); view(2); box on; end 
        if pstyle==3; p.pdeo.grid.plot(u,'LineStyle','none'); view(3);end        
        if pstyle==2; colorbar; end
      end     
      set(gca,'FontSize',p.plot.fs);
      xlabel('x'); ylabel('y');
      try colormap(p.plot.cm); catch, colormap cool; end;
      if ~strcmp(p.plot.axis,'tight'); axis(p.plot.axis); end 
      box on;
  case 3 %%%%%%%%%%%%%%%%%%%%%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    if isfield(p.plot,'ng')==0; p.plot.ng=20; end
    if pstyle==2; is=0; isn=0; % become 1 if lev or levn are passed via strings
       noa=nargin;
% adding extra options which are passed  via strings 
       for k=1:noa
         if ischar(varargin{k})==1
           switch lower(varargin{k})
             case 'lev'; p.plot.lev=varargin{k+1};is=1;
             case 'levc'; p.plot.levc=varargin{k+1};   
             case 'alpha'; p.plot.alpha=varargin{k+1};
             case 'ng'; p.plot.ng=varargin{k+1};
             case 'levn'; p.plot.levn=varargin{k+1};isn=1;
           end
         end
       end
% if field does not ex. set to 0.
       if isfield(p.plot,'lev')==0; p.plot.lev=0; end
       if isfield(p.plot,'levn')==0; p.plot.levn=0; end
       if isfield(p.plot,'levc')==0; p.plot.levc=0; end
% if lev is passed  via string ignore levn and vice versa.
       if is==1; p.plot.levn=0; end; if isn==1; p.plot.lev=0; end
% if lev and levn are both given, ignore levn.
       if iszero(p.plot.levn)==0 && iszero(p.plot.lev)==0; p.plot.levn=0; end
% if neither lev nor levn are given, set levn=1/4.
       if iszero(p.plot.lev)==1 && iszero(p.plot.levn)==1; p.plot.levn=1/4; end
% set levels, when levn is used
       if iszero(p.plot.levn)==0; k=p.plot.levn; 
          if k>=1; p.plot.lev=zeros(1,k); 
             for mm=1:k; p.plot.lev(mm)=min(u)+(max(u)-min(u))/(k+1)*mm; end; 
          else
             if k>0
               p.plot.lev=zeros(1,2); p.plot.lev(1)=min(u)+(max(u)-min(u))*k; p.plot.lev(2)=min(u)+(max(u)-min(u))*(1-k); 
             else
               k=-k; p.plot.lev=min(u)+(max(u)-min(u))*k;     
             end
           end
        end
% if levc does not ex. set default and if exist and is given via strings
% like 'r', convert to rgb code
        if iszero(p.plot.levc); p.plot.levc=[0 0 1; 1 0 0];
        else 
          if iscell(p.plot.levc); cell3d=p.plot.levc; 
             l=length(p.plot.levc); p.plot.levc=zeros(l,3); 
             for i=1:l; p.plot.levc(i,:)=pde2pathcolors(cell3d{i}); end
          else
             if ischar(p.plot.levc); p.plot.levc=pde2pathcolors(p.plot.levc); end
          end           
        end
% if two colors are given in levc, do color movement
        if size(p.plot.levc,1)==2 || size(p.plot.levc,2)==2
            nc=length(p.plot.lev); 
            if nc==1; p.plot.levc=p.plot.levc(1,:); 
            else
              c1=p.plot.levc(1,:); c2=p.plot.levc(2,:); d1=c2(1)-c1(1); d2=c2(2)-c1(2); d3=c2(3)-c1(3); 
              cm3d=zeros(nc,3); 
              for i=0:nc-1; cm3d(i+1,1)=c1(1)+d1*i/(nc-1); cm3d(i+1,2)=c1(2)+d2*i/(nc-1); cm3d(i+1,3)=c1(3)+d3*i/(nc-1); end
              p.plot.levc=abs(cm3d); 
            end
        end
        if isfield(p.plot,'alpha')==0; p.plot.alpha=1; end
        figure(wnr); clf; isoplot(p,u(n0:n1));  %plot solution
% convert colors from string to rgb code (used for creating a title)
if ischar(p.plot.levc); p.plot.levc=pde2pathcolors(p.plot.levc); end
if iscell(p.plot.levc); l=[]; for i=1:length(p.plot.levc); l=[l;  pde2pathcolors(p.plot.levc{i})]; end; p.plot.levc=l; end
%creating a title
l=length(p.plot.lev); 
if l==1; t=['level=\color[rgb]{',num2str(p.plot.levc(1,:)),'}',num2str(p.plot.lev,3)]; 
else
    t=['levels=\color[rgb]{',num2str(p.plot.levc(1,:)),'}',num2str(p.plot.lev(1),3)]; 
    for i=2:l
        t=[t,['\color{black}, \color[rgb]{',num2str(p.plot.levc(i,:)),'}',num2str(p.plot.lev(i),3)]];         
    end
end
     end % end of pstyle==2
if pstyle==1; slplot(p,p.plot.ng,u,p.plot.fs);  % p.pdeo.grid.plotSlices(u); 
    try colormap(p.plot.cm); catch,colormap cool;  end;  
    axis image;  colorbar('southoutside');  
end
if pstyle==3; faceplot(p,u); %p.pdeo.grid.plotFaces(u,'LineStyle','none');
    try colormap(p.plot.cm); catch, colormap cool; end;  
    colorbar('southoutside'); 
end %set(gca,'FontSize',p.plot.fs);
xlabel(''); ylabel(''); zlabel(''); 
end

if isfield(p.plot,'fs')==0; p.plot.fs=16; end
if p.nc.neq==1; t1=['u at ' tname]; % if we have a one component system
else t1=['u',num2str(cnr),' at ' tname]; 
end
% write title
if size(po,1)~=1 % title writing for 1D solutions is done in the 1 D part above  
if size(po,1)==3 && pstyle==2
    title({t1,t},'fontsize',p.plot.fs); % write point name and levels
else
    title(t1); % write point name
end
end
set(gca,'fontsize',p.plot.fs);
end % function plotsol ending

function anz=iszero(arg)% internal function, check if a structur field is zero (anz=1) or not (anz=0)
anz=0;
if isa(arg,'double');if max(size(arg))==1; if arg==0;anz=1; end; end; end
end