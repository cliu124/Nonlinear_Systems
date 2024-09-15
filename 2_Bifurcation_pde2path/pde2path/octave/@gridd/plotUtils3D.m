classdef plotUtils3D < handle
    
    methods(Access=public)        
        function plotFaces(obj,f,varargin)
            % triangle 
            indx=find(obj.e(4,:)==0); 
            switch nargin
                case 1
                    c =[0.9 0.9 0.9];
                otherwise 
                    if isempty(f)
                         c =[0.9 0.9 0.9];
                    else
                        c=reshape(f(obj.e(1:3,indx)),3,length(indx));
                    end                
            end
            x=reshape(obj.p(1,obj.e(1:3,indx)),3,length(indx));            
            y=reshape(obj.p(2,obj.e(1:3,indx)),3,length(indx));            
            z=reshape(obj.p(3,obj.e(1:3,indx)),3,length(indx));  
          %  patch(x,y,z,c,varargin{:},'EdgeColor','none');
          patch(x,y,z,c,varargin{:}); 
             
            % squares
            
            indx=find(obj.e(4,:)>0);  
            if isempty(indx)
                % there are NO sqares in the boundary, it must be a
                % tetraeder mesh!
                view(3)
                axis equal
                if nargin>2 && ~isempty(f) 
                    colorbar
                end
                return
            end
            switch nargin
                case 1
                    c =[0.9 0.9 0.9];
                otherwise
                    if isempty(f)
                         c =[0.9 0.9 0.9];
                    else
                        c=reshape(f(obj.e(1:4,indx)),4,length(indx));     
                    end
            end
            x=reshape(obj.p(1,obj.e(1:4,indx)),4,length(indx));            
            y=reshape(obj.p(2,obj.e(1:4,indx)),4,length(indx));            
            z=reshape(obj.p(3,obj.e(1:4,indx)),4,length(indx)); 
             
            patch(x,y,z,c,varargin{:});                  
            view(3)
            axis equal
            if nargin>2 && ~isempty(f) 
                colorbar
            end
        end
        
        function plottetras(obj,varargin)
            % triangle 
            nt=size(obj.t,2);  c =[0.9 0.9 0.9];
            indx=1:nt; 
            x=reshape(obj.p(1,obj.t(1:3,indx)),3,length(indx));            
            y=reshape(obj.p(2,obj.t(1:3,indx)),3,length(indx));            
            z=reshape(obj.p(3,obj.t(1:3,indx)),3,length(indx));  
            patch(x,y,z,c,varargin{:}); hold on; 
            x=reshape(obj.p(1,obj.t(2:4,indx)),3,length(indx));            
            y=reshape(obj.p(2,obj.t(2:4,indx)),3,length(indx));            
            z=reshape(obj.p(3,obj.t(2:4,indx)),3,length(indx));  
            patch(x,y,z,c,varargin{:}); 
            x=reshape(obj.p(1,obj.t([1 2 4],indx)),3,length(indx));            
            y=reshape(obj.p(2,obj.t([1 2 4],indx)),3,length(indx));            
            z=reshape(obj.p(3,obj.t([1 2 4],indx)),3,length(indx));  
            patch(x,y,z,c,varargin{:}); 
            view(3);             
        end
        
        
        function plotIso(obj,u,lev,varargin) 
        % iso-surfaceplot 
        % obj.plotIso(y,iso,options)
        % Plots the iso-surface of the data y to the level iso
        % Optional arguments are NGP (Number of Grid Points of the iso
        % surface) and the FontSize of the plot.
        % This method bases on the request and the prototype code
        % of Hannes Uecker. Thanks!
        % (C) 2014 by Uwe Prüfert
        
        
            ng=20;
            fs=14;
            color=[0.6 0.3 0.8
                     0.8 0.8 0.8
                     0.5 0.5 1
                     0.6 0.6 1
                     0.2 0.7 0.8];
            
            if all(lev>max(u)) && all(lev < min(u))
                fprintf('Empty level set\n');
            end
       
            gp=obj.p; 
            
            x1=min(gp(1,:)); 
            x2=max(gp(1,:)); 
            y1=min(gp(2,:)); 
            y2=max(gp(2,:)); 
            z1=min(gp(3,:));
            z2=max(gp(3,:)); 
            
            xv=linspace(x1,x2,ng); 
            yv=linspace(y1,y2,ng); 
            zv=linspace(z1,z2,ng); 
            
            [X,Y,Z]=meshgrid(xv,yv,zv);
            
            up=obj.p3interpol(X,Y,Z,u,gp(1,:),gp(2,:),gp(3,:));
            
            if length(lev)>5
                warning('PLOTUTIL3D:TOOMANYLEVELS',...
                    'The number of iso surfaces to plot is restricted to five.')
            end
            hold on
            for k=1:min(length(lev),5)                
                ip=patch(isosurface(X,Y,Z,up,lev(k)));               
                isonormals(X,Y,Z,up,ip);
                set(ip,'FaceColor',color(k,:),'EdgeColor','none');
            end
                
            view(3); axis([x1 x2 y1 y2 z1 z2]);
            camlight;
            lighting phong;
            xlabel('x','FontSize',fs); 
            ylabel('y','FontSize',fs); 
            zlabel('z','FontSize',fs); 
            grid on;
            title(['Level=' mat2str(lev,3)],'FontSize',fs);
            set(gca,'FontSize',fs); 
            plot3(x1,y1,z1);
            plot3(x2,y2,z2);
            axis equal
            hold off
            if nargin>3
                % additional options                
                set(ip,varargin{:});
            end
        end
        
        function cutawayPlot(obj,varargin)
        % cutawayPlot - cuts away all points x<X, y<Y z<Z and
        % plots the object.
        % obj.cutawayPlot(X,Y,Z) 
        % Two of the arguments X,Y,Z may be empty.
        % obj.cutawayPlot(X,Y,Z,y) plots the object with the data y
        % (c) 2014 by Uwe Prüfert
            
            
            % make a copy of all data
          
            p=obj.p;
            
            % we do not need the edge number information here...
            e=obj.e(1:4,:);
            t=obj.t;
            
            % find all point with x<=X etc...
            
            switch length(varargin)
                case 3
                    X=varargin{1}; Y=varargin{2};  Z=varargin{3}; 
                    X2=[]; Y2=[];  Z2=[]; y=[];                    
                case 4
                    X=varargin{1};
                    Y=varargin{2};
                    Z=varargin{3};
                    X2=[];
                    Y2=[];
                    Z2=[];
                    y=varargin{4};
                case 6
                    X=varargin{1};                    
                    Y=varargin{2};
                    Z=varargin{3};
                    X2=varargin{4};
                    Y2=varargin{5};
                    Z2=varargin{6};                   
                case 7
                   X=varargin{1};
                   Y=varargin{2};
                   Z=varargin{3};
                   X2=varargin{4};
                   Y2=varargin{5};
                   Z2=varargin{6};
                   y=varargin{7};                   
                otherwise
                    obj.wrongNumberInputs.throw
            end 
            if ~isempty(X)             
                indX=find(obj.p(1,:)<X);
            else
                indX=[];
            end
            if ~isempty(Y)             
                indY=find(obj.p(2,:)<Y);
            else
                indY=[];
            end
            if ~isempty(Z)             
                indZ=find(obj.p(3,:)<Z);
            else
                indZ=[];
            end
            if ~isempty(X2)             
                indX2=find(obj.p(1,:)>X2);
            else
                indX2=[];
            end
            if ~isempty(Y2)             
                indY2=find(obj.p(2,:)>Y2);
            else
                indY2=[];
            end
            if ~isempty(Z2)             
                indZ2=find(obj.p(3,:)>Z2);
            else
                indZ2=[];
            end
            indX=unique([indX indY indZ indX2 indY2 indZ2]);
            if isempty(indX)  
                MException('grid3Dpr:cutawayPlot:EmptySet',...
                    'The cut is not well defined, check the arguments.').throw
                
            end
            
            % remove points
            numberOfRemovedPoints=length(indX);
            fprintf([num2str( numberOfRemovedPoints),' points are in the cutaway.\n'])
 
            for k=indX
                % look for removed points in elements and set the
                % position of removed points to zero.
                % indX(k) found in r-th row and c-th column
                 
                [r,c]=find(e==k);  
                
                for k2=1:length(r)
                    e(r(k2),c(k2))=0;
                end
                [r,c]=find(t==k);  
                
                for k2=1:length(r)
                    t(r(k2),c(k2))=0;
                end                
            end
            
            % now all zeros t can be removed and 
            % all t with less than 3 non-zeros
            % all three and four non-zeros become new edges
            k=1;
            % remove deformed edges
            while true
                if k>size(e,2)
                    break;
                end
                switch length(find(e(:,k)>0))
                    case {3,4}
                        % edge element triangle or                      
                        % edge element rectangle 
                        % stays in e so we skip the removing                        
                        k=k+1;
                    otherwise % 0 1 2
                        % to removed: patches with removed points
                        % do not increase k because k + 1 "jumps" to k
                        % "automaticly"
                        e(:,k)=[];
                end
            end            
            
            % add former elements (now with three or four points insted 
            % of six to edges
            for k=1:obj.nElements                
                switch length(find(t(:,k)>0))
                    case 3
                        % edge element triangle or                      
                        % edge element rectangle 
                        % stays in e so we skip the removing 
                        e=[e ,[t((t(:,k)>0),k);0]]; 
                    case 4
                        % edge element triangle or                      
                        % edge element rectangle 
                        % stays in e so we skip the removing      
                        id=find(t(:,k)>0); 
                        % correct the order of points in the "cut in x or
                        % cut in y" case
                        e=[e ,t(id([1 2 4 3]),k)];                         
                    otherwise % 0 1 2
                        %  do nothing
                end
            end
 
            % create a grid3Dpr object and fill it with the OLD points and
            % NEW edges. The t is a dummy, becaus grid3Dpr.display needs a
            % nonempty t. It 's only to prevent "mistyrious errors" 
            g=grid3Dpr;
            g.setPET(p,e,1);
                       
 
            % use the plotFaces
            if isempty(y)
                g.plotFaces(); 
            else
                g.plotFaces(y);
            end
        end
        
    function cutplot(obj,y,cut) % HU 
    % cuts away all points x<X (and x>X2), y<Y, z<Z and plots u on remainder
    % obj.cutawayPlot(u,cut), where cut=[X,Y,Z] or cut=[X,Y,Z,X2Y2,Z2] 
    % mod of Uwe's cutawayplot 
        p=obj.p;         

    % we do not need the edge number information here...
    e=obj.e(1:4,:); t=obj.t;
    switch length(cut) % find all point with x<=X etc...
        case 3;  X=cut(1); Y=cut(2);  Z=cut(3); X2=[]; Y2=[];  Z2=[]; 
        case 6; X=cut(1); Y=cut(2); Z=cut(3); X2=cut(4); Y2=cut(5); Z2=cut(6); 
    end 
    if ~isempty(X); indX=find(obj.p(1,:)<X); else indX=[]; end
    if ~isempty(Y); indY=find(obj.p(2,:)<Y); else indY=[]; end
    if ~isempty(Z); indZ=find(obj.p(3,:)<Z); else indZ=[]; end
    if ~isempty(X2); indX2=find(obj.p(1,:)>X2); else indX2=[]; end
    if ~isempty(Y2); indY2=find(obj.p(2,:)>Y2); else indY2=[]; end
    if ~isempty(Z2); indZ2=find(obj.p(3,:)>Z2); else indZ2=[]; end
    indX=unique([indX indY indZ indX2 indY2 indZ2]);       
    % remove points
   % numberOfRemovedPoints=length(indX); fprintf([num2str( numberOfRemovedPoints),' points in the cutaway.\n'])

    for k=indX % look for removed points in elements and set the  
% position of removed points to zero. indX(k) found in r-th row and c-th column
      [r,c]=find(e==k);  
      for k2=1:length(r); e(r(k2),c(k2))=0; end
      [r,c]=find(t==k);  
      for k2=1:length(r); t(r(k2),c(k2))=0; end                
    end
% now all zeros t can be removed and all t with less than 3 non-zeros
% all three and four non-zeros become new edges
    k=1;
    % remove deformed edges
    while true
        if k>size(e,2); break; end
        switch length(find(e(:,k)>0))
            case {3,4}; % edge element triangle or edge element rectangle stay 
                k=k+1;
            otherwise % 0 1 2; to removed: patches with removed points
                % do not increase k because k + 1 "jumps" to k "automaticly"
                e(:,k)=[];
        end
    end            

    % add former elements (now with three or four points insted of six to edges
    for k=1:obj.nElements                
        switch length(find(t(:,k)>0))
            case 3; % edge element triangle or edge element rectangle 
                % stays in e so we skip the removing 
                e=[e ,[t((t(:,k)>0),k);0]]; 
            case 4;  
                id=find(t(:,k)>0); 
                % correct the order of points in the "cut in x or cut in y" case
                e=[e ,t(id([1 2 4 3]),k)];                         
            otherwise % 0 1 2
                %  do nothing
        end
    end
% create a grid3Dpr object and fill it with the OLD points and NEW edges. The t is a dummy, 
% becaus grid3Dpr.display needs a nonempty t. 
    g=grid3Dpr; g.setPET(p,e,1);    
% use the plotFaces
    if isempty(y); g.plotFaces(); 
    else; g.plotFaces(y);
    end
end % cutplot, HU 
        
        function un=getSlice(obj,u,slice,ng,fs)
            % getSlices gives back points and values of slices wrt. 3D
            % objects. Essentially the same as plotSlices
            x1=min(obj.p(1,:)); x2=max(obj.p(1,:));
            y1=min(obj.p(2,:)); y2=max(obj.p(2,:)); 
            z1=min(obj.p(3,:)); z2=max(obj.p(3,:)); 
            n=10;
            x=linspace(0,1,n);
            [X,Y]=meshgrid(x,x);
            X=[reshape(X,n^2,1) reshape(Y,n^2,1) 0.5*ones(n^2,1)];
            un=obj.p3interpol(X(:,1),X(:,2),X(:,3),u,obj.p(1,:),obj.p(2,:),obj.p(3,:));
            size(un) 
            un=reshape(un,n,n);
        end
        
        
        
        function plotSlices(obj,u,xslice,yslice,zslice,varargin)
        % plotSlices - plots slices of a 3D function over 3D grid object.    
        % plotSlices(obj,u)  
        % plotSlices(obj,u,xslice,yslice,zslice)  
        % plotSlices(obj,u,xslice,yslice,zslice[,options])    
        % xslice,yslice,zslice are vectors given the coordinates of
        % the slices.
        % Every of the arguments xslice,yslice,zslice can be empty. 
        % In this case it will be set to a mean value. 
        % This method bases on the request and the prototype code
        % of Hannes Uecker. Thanks!
        % (C) 2014 by Uwe Prüfert   
        
            x1=min(obj.p(1,:)); x2=max(obj.p(1,:));
            y1=min(obj.p(2,:)); y2=max(obj.p(2,:)); 
            z1=min(obj.p(3,:)); z2=max(obj.p(3,:)); 
            
            % 
            ng=40;
            fs=14;
            
            switch nargin
                case {1,2}
                    %only u given.
                    xslice=0.5*(x1+x2); 
                    yslice=0.5*(y1+y2);
                    zslice=0.5*(z1+z2);
                otherwise 
                    % full set of parameters given       
            end
            
            if isempty(xslice)&&isempty(yslice)&&isempty(zslice)
                xslice=0.5*(x1+x2);             
                yslice=0.5*(y1+y2);            
                zslice=0.5*(z1+z2);
            end                  
            xv=linspace(x1,x2,ng); yv=linspace(y1,y2,ng); zv=linspace(z1,z2,ng); 
            
            xv=unique([xv , xslice]);
            yv=unique([yv , yslice]);
            zv=unique([zv , zslice]);
            
            [X,Y,Z]=meshgrid(xv,yv,zv);
            up=obj.p3interpol(X,Y,Z,u,obj.p(1,:),obj.p(2,:),obj.p(3,:));                   
            
            sl=slice(X,Y,Z,up,xslice,yslice,zslice);
            view(3); axis([x1 x2 y1 y2 z1 z2]); 
            xlabel('x','fontsize',fs); ylabel('y','fontsize',fs); 
            zlabel('z','fontsize',fs);set(gca,'FontSize',fs) 
            set(sl,'EdgeColor','none','FaceColor','Interp'); 
            if ~isempty(varargin)
                set(sl,varargin{:});
            end
            colorbar
            hold on
            plot3(x1,y1,z1);
            plot3(x2,y2,z2);
            axis equal
        end 
        
        function moveMesh(obj,xshift,yshift,zshift)
            % moveMesh 
            % shifts all points in the mesh with xshift,yshift zshift.
            % Note, that only the points are changed.  
            % (c) 2013 by Uwe Prüfert
            obj.p(1,:)=obj.p(1,:) + xshift;
            obj.p(2,:)=obj.p(2,:) + yshift;
            obj.p(3,:)=obj.p(3,:) + zshift;
        end
        
        function identifyBoundarySegment(obj,bdno)
            % identifyBoundary - plots the boundary segments in
            % colors.
            % obj.identifyBoundary every segment has its own color,
            % starting from blue over green to red.
            % obj.identifyBoundary(BDNR) plots the BDNR-th boundary
            % segment in red, all other boundary segments in blue.
            % (c) 2014 by Uwe Prüfert
            if nargin == 2
                c= zeros(1,obj.nEdges);
                c(obj.e(5,:)==bdno)=1;
            else
                n=ceil(obj.nBoundarySegments/3);
                if n >3
                    warndlg({['Large number of boundary segments',...
                        ' cannot display it this way.'];...
                        'Use identifyBoundarySegment(noSegment) insted'},...
                        'Warning')
                    return
                end
                for k=1:obj.nBoundarySegments
                     subplot(n,3,k)
                     title(['Boundary segment no ', num2str(k),])
                     obj.identifyBoundarySegment(k)
                     axis off
                end
                return
            end
            indx=find(obj.e(4,:)==0); 
            x=reshape(obj.p(1,obj.e(1:3,indx)),3,length(indx));            
            y=reshape(obj.p(2,obj.e(1:3,indx)),3,length(indx));            
            z=reshape(obj.p(3,obj.e(1:3,indx)),3,length(indx));  
             
            patch(x,y,z,c(indx));             
            % squares
            
            indx=find(obj.e(4,:)>0);  
            x=reshape(obj.p(1,obj.e(1:4,indx)),4,length(indx));         
            y=reshape(obj.p(2,obj.e(1:4,indx)),4,length(indx));            
            z=reshape(obj.p(3,obj.e(1:4,indx)),4,length(indx)); 
            
             
            patch(x,y,z,c(indx));     
            % patch(x,y,z,'r');
            view(3)
            axis equal 
        end
        
        function indx=pointToElementIndex(obj,pt)
            % pointToTriangleIndex-  gives back the index of the triangle
            % containing the point pt pt must be a vector of length
            % two.
            % indx=isPointInTriangle(obj,pt)
           
            b= obj.isPointInElement(pt);
            indx=find(b);            
        end
        
        function b=isPointInDomain(obj,pt)
            % isPointInDomain - true if the point pt is in at least
            % one triangle of the triangulation
            % 
            b=any(obj.isPointInElement(pt));
        end
    end
    
    methods(Static,Access=protected,Hidden)    
        function un=p3interpol(xn,yn,zn,u,x,y,z)
            % interpolate u def. on x,y,z to new mesh given in xn,yn,zn 
            xv=reshape(x, size(x,1)*size(x,2), 1);
            yv=reshape(y, size(y,1)*size(y,2), 1);
            zv=reshape(z, size(z,1)*size(z,2), 1);
            uv=reshape(u, size(u,1)*size(u,2), 1);
         %   F=scatteredInterpolant(xv,yv,zv,uv);    un=F(xn,yn,zn);
         un= griddata3 (xv, yv, zv, uv, xn, yn, zn)
        end
    end
end