function [h,XXr,YYr]=mypdeplot(varargin)
% MYPDEPLOT: some code taken from pdetoolbox pdeplot to plot and RETURN contours 
% (to extract spiral tip) 
%
%   See also: PDEPLOT3D, pde.PDEMODEL, PDECONT, PDEGPLOT, pde.FEMesh, PDESURF
%       Copyright 1994-2016 The MathWorks, Inc.     
% Default values
xydata=[];  zdata=[]; zstyle='continuous'; 
cmap='cool'; intool='off'; xygrid='off'; mh='off';
title=''; levels=10; tn=[]; a2=[]; a3=[]; znodedata=0; 

% Recover Param/Value pairs from argument list
nargs=nargin;  p=varargin{1}; e=varargin{2}; t=varargin{3}; varargin={varargin{4:end}};
numVarArgs=nargs - 3; xystyle='off'; cont='on'; plotxy=1; plotz=0; intool=0; 
for ii=1:2:numVarArgs
  Param=varargin{ii}; Value=varargin{ii+1};
  if ~ischar(Param)
    error(message('pde:pdeplot:ParamNotString'))
  elseif size(Param,1)~=1
    error(message('pde:pdeplot:ParamNumRowsOrEmpty'))
  end
  switch lower(Param)
  case 'xydata'
    xydata=Value;
  case 'levels'
    levels=Value;
    if isempty(levels); levels=10; end
  end
end

ntri=size(t,2); nnode=size(p,2);
if plotxy
  xys=size(xydata);
  if xys(1)==nnode; xynodedata=1; xydata=xydata(:,1);
  elseif xys(2)==ntri; xynodedata=0; xydata=xydata(1,:);
  elseif xys(2)==nnode; xydata=xydata'; xynodedata=1; xydata=xydata(:,1);
  elseif xys(1)==ntri; xydata=xydata'; xynodedata=0; xydata=xydata(1,:);
  else; error(message('pde:pdeplot:xydataLength'))
  end
end

ax=newplot; hold on; view(2); colormap(cmap); hh=[]; 

% case: contour plot
if strcmp(cont,'on') && plotxy
  zdata=xydata;
  xymin=min(min(xydata)); xymax=max(max(xydata)); 
  zmin=min(min(zdata)); zmax=max(max(zdata));
  if xymax==xymin, xymax=xymin+1; end
  if zmax==zmin, zmax=zmin+1; end
  if numel(levels)==1
    n=levels; lmin=(n*xymin+xymax)/(n+1); lmax=(xymin+n*xymax)/(n+1);
    levmin=xymin; levmax=xymax;  zlmin=(n*zmin+zmax)/(n+1);
    zlmax=(zmin+n*zmax)/(n+1); lev=linspace(lmin,lmax,n);  zlev=linspace(zlmin,zlmax,n);
  else
    levels=sort(levels); n=length(levels); lmin=levels(1);    lmax=levels(n);
    zlmin=lmin;  zlman=lmax; lev=levels; zlev=levels; levmin=xymin; levmax=xymax;
  end

  cm=colormap;  ncm=size(cm,1);
  icm=floor(((lev-levmin)/(levmax-levmin))*(ncm-1)+0.5)+1;
  if max(icm)>ncm || min(icm)<1
    icmindx=find(icm<=ncm & icm>=1);
    icm=icm(icmindx);
  end
  ccm=cm(icm,:);

  % Ensure that overlayed contour is drawn on top of surface plot
  if ~strcmp(xystyle,'off') && ~plotz;  set(gca,'SortMethod','childorder'); end

  if strcmp(xygrid,'on')
  else  
     if size(t, 1)==4;  tlq=t;
     else;  tlq=[t([1, 4, 6], :), t([4, 2, 5], :), t([5, 3, 6], :), t([4, 5, 6], :)];
     end      
    nt=size(tlq,2); zt=reshape(zdata(tlq(1:3,:)),3,nt);
    xyt=reshape(xydata(tlq(1:3,:)),3,nt); ztmax=max(zt); ztmin=min(zt); XX=[]; YY=[]; ZZ=[];
    for j=1:length(lev)
      jlev=zlev(j);
      it=find(ztmin<=jlev & ztmax>=jlev);
      if size(it)
        z1=zt(1,it); z2=zt(2,it); z3=zt(3,it);
        a21=zeros(1,length(it));
        itt=find(z2~=z3);       % This kludge is to avoid the warning message
        a21(itt)=(jlev-z3(itt))./(z2(itt)-z3(itt));
        itt=find(z2==z3);
        a21(itt)=NaN*ones(size(itt));
        a32=zeros(1,length(it));
        itt=find(z3~=z1);
        a32(itt)=(jlev-z1(itt))./(z3(itt)-z1(itt));
        itt=find(z3==z1);
        a32(itt)=NaN*ones(size(itt));
        a13=zeros(1,length(it));
        itt=find(z1~=z2);
        a13(itt)=(jlev-z2(itt))./(z1(itt)-z2(itt));
        itt=find(z1==z2);
        a13(itt)=NaN*ones(size(itt));

        a2=NaN*ones(2,length(it));  a3=NaN*ones(2,length(it));
        ii=ones(1,length(it));          % 1+the number of points found so far

        itt=find(a21>=0 & a21<=1);      % On side 1
        a2(ii(itt)+2*(itt-1))=a21(itt);
        a3(ii(itt)+2*(itt-1))=1-a21(itt);
        ii(itt)=ii(itt)+ones(size(itt));
        itt=find(a32>=0 & a32<=1);      % On side 2
        a2(ii(itt)+2*(itt-1))=zeros(size(itt));
        a3(ii(itt)+2*(itt-1))=a32(itt);
        %  ii(itt)=ii(itt)+ones(size(itt));
        itt=find(a13>=0 & a13<=1);      % On side 3
        % This must be the second endpoint
        a2(2,itt)=1-a13(itt);
        a3(2,itt)=zeros(size(itt));

        X=[(1-a2(1,:)-a3(1,:)).*p(1,tlq(1,it))+ ...
                a2(1,:).*p(1,tlq(2,it))+a3(1,:).*p(1,tlq(3,it)); ...
            (1-a2(2,:)-a3(2,:)).*p(1,tlq(1,it))+ ...
                a2(2,:).*p(1,tlq(2,it))+a3(2,:).*p(1,tlq(3,it)); ...
            NaN*ones(size(it))];
        Y=[(1-a2(1,:)-a3(1,:)).*p(2,tlq(1,it))+ ...
                a2(1,:).*p(2,tlq(2,it))+a3(1,:).*p(2,tlq(3,it)); ...
            (1-a2(2,:)-a3(2,:)).*p(2,tlq(1,it))+ ...
                a2(2,:).*p(2,tlq(2,it))+a3(2,:).*p(2,tlq(3,it)); ...
            NaN*ones(size(it))];
        Z=[jlev*ones(size(it)); jlev*ones(size(it)); NaN*ones(size(it))];
        X=X(:);
        Y=Y(:);
        Z=Z(:);

        nxx=size(XX,1);
        nx=size(X,1);
        if nxx>nx
          nn=NaN*ones(nxx-nx,1);
          X=[X; nn];
          Y=[Y; nn];
          Z=[Z; nn];
        elseif nxx<nx
          nn=NaN*ones(nx-nxx,size(XX,2));
          XX=[XX; nn];
          YY=[YY; nn];
          ZZ=[ZZ; nn];
        end
        XX=[XX, X];
        YY=[YY, Y];
        ZZ=[ZZ, Z];
      end                               % size(it)
    end

    % plot geometry boundaries:        
    hold on

    if ~plotz
      ZZ=zeros(size(XX));
    end
    for i=1:size(XX,2)
      if strcmp(xystyle,'off')
      % Colored contours:
        hndl(i)=line(XX(:,i),YY(:,i),ZZ(:,i),...
            'Parent',ax,...
            'color',ccm(i,:)); %HU
        XXr=XX(~isnan(XX)); YYr=YY(~isnan(YY)); % HU
      else
      % Overlayed black contours:
        contc='k';
        hndl(i)=line(XX(:,i),YY(:,i),zeros(size(XX(:,i))),...
            'Parent',ax,...
            'color',contc);
      end
    end
    h1=pdeplot(p,e,[]);
    hold off

    if(exist('hndl', 'var'))
      hhc=[hndl(:);h1(:)];
    else
      hhc=h1(:);
    end
  end
  hh=[hh; hhc];  
end % if strcmp(cont,'on') && plotxy,
% Finally, set the axes title
col='k';
set(get(ax,'Title'), 'Color',col, 'String',title, 'Visible','on')
hold off
h=hh; %HU
end