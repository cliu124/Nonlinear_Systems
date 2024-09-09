classdef grid2D < gridd    
 % properties(SetAccess=protected, GetAccess=public); nPointsInElements double=3;    end
 % properties(Constant=true); spaceDimension double=2; end
    methods(Access=public) 
        % constructor  and copy
        function obj=grid2D()
            obj.nPointsInElements=3; 
            switch nargin
                case 0
                    % empty object
                otherwise
                    obj.wrongNumberInputs.throw;
            end          
        end  
        
function pts2dom(obj,P) % point-cloud to generate PTE 
dt=delaunayTriangulation(P'); obj.p=dt.Points'; obj.e=dt.freeBoundary'; obj.t=dt.ConnectivityList';
%T=obj.t; T=T(1:4,:)'; X=obj.p'; figure(10); clf; tetramesh(T,X); pause 
obj.t(5,:)=ones(1,size(obj.t,2));  obj=setids1(obj); 
end   
        % basic geometries: squares l-shapes, circles                      
        function lshape(obj,hmax)
            % Lshapeq
            % obj.lshape([hmax])
             obj.p=[0 0.5 1 0   0.5 1   0 0.5;...
                     0 0   0 0.5 0.5 0.5 1 1 ];
                 
            obj.e=[1   2   3 6 5 8 7   4;
                    2   3   6 5 8 7 4   1;...
                     0   0.5 0 0 0 0 0   0.5;...
                     0.5 1   1 1 1 1 0.5 1;...
                     1   1   2 3 4 5 6   6;...
                     1   1   1 1 1 1 1   1;...
                     0   0   0 0 0 0 0   0];
                 
            obj.t=[1 2 2 3 4 5 ;...
                     2 5 3 6 5 8;...
                     4 4 5 5 7 7;...
                     1 1 1 1 1 1];
                 
            switch nargin
                case 1
                case 2
                    if ischar(hmax)               
                        hmax=str2double(hmax);
                    end
                    h=sqrt(2);
                    while h > hmax
                        obj.refineMesh;
                        h=h/2;
                    end                 
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end
        
        function unitSquare(obj,hmax)
           % a version of unitsquare that creates criss-cross
            % triangulations
            % obj.unitsquare[hmax])           
            % (c) 2013 Uwe Prüfert
            
            obj.p=[0 1 1 0;...
                   0 0 1 1];
                 
            obj.e=[1 2 3 4;
                     2 3 4 1;...
                     0 0 0 0;...
                     1 1 1 1;...
                     1 2 3 4;...
                     1 1 1 1;...
                     0 0 0 0];
                 
            obj.t=[1 2;...
                     2 3;...
                     4 4;...
                     1 1];
                 
            switch nargin
                case 1
                    % nothing to do hmax=Inf
                case 2
                    % ok
                    if ischar(hmax)               
                        hmax=str2double(hmax);
                    end
                    h=sqrt(2);
                    while h > hmax
                        obj.refineMesh;
                        h=h/2;
                    end
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end        
        
        
         function mySquare(obj,lx,ly,hmax)
           % oct version, don't use for bad aspect ratios
            obj.p=[-lx lx lx -lx;...
                   -ly -ly ly ly];   
                 
            obj.e=[1 2 3 4;
                   2 3 4 1;...
                   0 0 0 0;...
                   1 1 1 1;...
                   1 2 3 4;...
                   1 1 1 1;...
                   0 0 0 0];
                 
            obj.t=[1 2;...
                   2 3;...
                   4 4;...
                   1 1];
                 
             h=sqrt(lx/2);
             while h > hmax; obj.refineMesh;  h=h/2;
         
            end
        end        
        
        function doubleT(obj,scaleFactor)
            %doubleT
            % Creates a grid for a double-t shaped domain
            % doubleT()
            % doubleT(scale)
            % Called with no argument, the base hal lenght one and the
            % height is three. One can scale the whole domain, but not
            % the proportions.
            % There are 12 boundary segments.
            % To identify boundary segment number call
            % obj.identifyBoundarySegment 
            % (c) 2013 by Uwe Prüfert.
            
            switch nargin
                case 1
                    scale=1/3;
                case 2
                    % scale factor
                    if ~(isscalar(scaleFactor)&&isa(scaleFactor,'double'))
                        obj.wrongClass.throw;
                    end
                otherwise
                    obj.wrongNumberInputs.throw;
            end
            
            obj.p=[0 1 2 3 0 1 2 3 1 2 1 2 0 1 2 3 0 1 2 3;...
                     0 0 0 0 1 1 1 1 2 2 3 3 4 4 4 4 5 5 5 5];
             
            obj.e=[1   2   3   4 8 7   10  12  15 16  20  19  18  17 13 14  11   9   6  5;...
                     2   3   4   8 7 10  12  15  16 20  19  18  17  13 14 11  9    6   5  1;...
                     0   1/3 2/3 0 0 0   1/3 2/3 0  0   0   1/3 2/3 0  0  0   1/3  2/3 0  0;...
                     1/3 2/3 1   1 1 1/3 2/3 1   1  1   1/3 2/3 1   1  1  1/3 2/3  1   1  1;...
                     1   1   1   2 3 4   4   4   5  6   7   7   7   8  9  10  10   10  11 12];

            obj.t=[1 2 2 3 3 4 6 7  9  10 11 12 13 14 14 15 15 16;...
                     2 6 3 7 4 8 7 10 10 12 12 15 14 18 15 19 16 20;...
                     5 5 6 6 7 7 9 9  11 11 14 14 17 17 18 18 19 19];
            obj.p = obj.p/3; 
        end
        
        function unitCircle(obj,hmax)
            % obj.unitcircle
            % Uses drid2D.circle 
            % obj.unitcircle()
            % obj.unitcircle(hmax)
            % (c) 2013 Uwe Prüfert
             
            switch nargin
                case 1
                    obj.circle(1,0.2);      
                case 2
                    obj.circle(1,hmax);
                otherwise
                    obj.wrongNumberInputs.throw;
            end          
                  
        end    
        
        function freeGeometry(obj,varargin)
            % obj.freeFormedGeometry(outerConstraint,[innerconstrain(s)])
            %(c) 2013 by Uwe Prüfert
            
            if nargin==1
                % no constraint given?
                obj.wrongNumberInputs.throw;
            elseif nargin==2
                X=varargin{1};  [xx,xy]=size(X);
                if min(xx,xy)~=2
                    obj.wrongFormat.throw;
                end
                if xy==2
                    X=X';
                end
                % distance of boundary points
                ds= sqrt(sum([(X(:,2:end)-X(:,1:end-1)).^2 (X(:,1)-X(:,end)).^2]));

                if ds(end) < 1e-6
                    % we assume now that X ist a closed curve and
                    % we do not need the doubled point
                    X(:,end)=[];   ds(end)=[];
                end

                % update number of points

                np=length(ds);        h=mean(ds);
                nx=ceil((max(X(1,:))-min(X(1,:)))/h)+1;
                ny=ceil((max(X(2,:))-min(X(2,:)))/h)+1;


                [xm,ym]=meshgrid(linspace(min(X(1,:))+h,max(X(1,:))-h,nx-2),...
                    linspace(min(X(2,:))+h,max(X(2,:))-h,ny-2));
                [nx,ny]=size(xm);
                xm=reshape(xm,1,nx*ny);     ym=reshape(ym,1,nx*ny);

                % the constraint is the index of boundary elements
                C= [1:np;2:np,1]';

                X=[X';[xm;ym]'];                    
               
                     dt=delaunayTriangulation(X,C);
                     p=dt.Points';  
                     t=dt.ConnectivityList';
                     io=isinside(1:size(t,2)); % HU to do

                % Detect all triangles inside the object
                % all points, including out side points
                % Select all points from inside triangles
                p=p(:,unique(t(:,io))');
                t=t(:,io);

                % New triangulation only with inside points
                % We need this iteration because of the renumbering
              
                     dt=delaunayTriangulation(p',C);
                     t=dt.ConnectivityList';
                     io=1:size(t,2); 
%                   
                obj.p=p;                    
                obj.e=C';

                obj.t=t(:,io);
                obj.t(4,:)=ones(1,size(obj.t,2)); 

                % the parameter s runs from 0 to 1 on the boundary
                sds=sum(ds);

                obj.e(3,1)=0;
                obj.e(4,1)=ds(1)/sds;
                obj.e(5,1)=1;

                % parametrize the boundary
                for k=2:np
                    obj.e(3,k)=obj.e(4,k-1);   
                    obj.e(4,k)=obj.e(4,k-1) + ds(k) /sds;
                    obj.e(5,k)=1; % bd segment number
                end
            else
                % a number of (interior) boundaries.
                % prepare data
                % outa "constraint"
                % Collect all boundary segments

                % We check the "polynomiality"

                dsmean=zeros(1,length(varargin));
                for k=1:length(varargin)
                    X=varargin{k};
                    [xx,xy]=size(X);
                    if min(xx,xy)~=2
                        obj.wrongFormat.throw;
                    end
                    if xy==2
                        % it should be a matrix x;y coordinate
                        X=X';
                    end
                    ds= sqrt(sum([(X(:,2:end)-X(:,1:end-1)).^2 (X(:,1)-X(:,end)).^2]));
                    dsmean(k)=mean(ds);
                end

                dsmin=min(dsmean);

                XX=[];
                XC=[];
                NBDS=[]; % the 5-th row in obj.e, 
                           % the number of boundary segment 
                XDS=[];  
                XE=[];
                for k=1:length(varargin)
                    X=varargin{k};
                    [xx,xy]=size(X);
                    if min(xx,xy)~=2
                        obj.wrongFormat.throw;
                    end
                    if xy==2
                        % it should be a matrix x;y coordinate
                        X=X';
                    end
                    ds= sqrt(sum([(X(:,2:end)-X(:,1:end-1)).^2 (X(:,1)-X(:,end)).^2]));
                    % Here prevent "huge edges"


                    if dsmean(k)/dsmin>3
                        XN=[];
                        % refine the polygonial
                        xp=1:dsmin*0.5:2;

                        for k1=2:size(X,2)
                            xtemp=interp1([X(1,k1-1) X(1,k1)],xp);
                            ytemp=interp1([X(2,k1-1) X(2,k1)],xp);
                            Xtemp=[xtemp;ytemp];
                            Xtemp(:,end)=[];
                            XN=[XN Xtemp];
                        end
                        xtemp=interp1([X(1,k1) X(1,1)],xp);
                        ytemp=interp1([X(2,k1) X(2,1)],xp);
                        Xtemp=[xtemp;ytemp];
                        Xtemp(:,end)=[];
                        XN=[XN Xtemp];
                        X=XN;
                    end                     
                    ds= sqrt(sum([(X(:,2:end)-X(:,1:end-1)).^2 (X(:,1)-X(:,end)).^2]));                     

                    if ds(end) < 1e-6
                        % we assume now that X ist a closed curve and
                        % we do not need the doubled point
                        X(:,end)=[];
                        ds(end)=[];
                    end
                    np=length(ds);

                    if k==1
                        C=[1:np;...
                            2:np,1];
                    else
                        C= [XC(1,end)+1:XC(1,end)+np;...
                              XC(1,end)+2:XC(1,end)+np,XC(1,end)+1];
                    end
                    NBDS=[NBDS,k*ones(1,np)];
                    XC=[XC C];
                    XX=[XX X];
                    XDS=[XDS ds];

                    sds=sum(ds);

                    e3=0;
                    e4=ds(1)/sds;                      

                    for k1=2:np                             
                        e3(k1)=e4(k1-1);   
                        e4(k1)=e4(k1-1) + ds(k1) /sds;                                                       
                    end                        
                    XE=[XE [e3;e4]];
                end


                h=mean(XDS);
                nx=ceil((max(XX(1,:))-min(XX(1,:)))/h)+1;
                ny=ceil((max(XX(2,:))-min(XX(2,:)))/h)+1;


                [xm,ym]=meshgrid(linspace(min(XX(1,:))+h,max(XX(1,:))-h,nx-2),...
                    linspace(min(XX(2,:))+h,max(XX(2,:))-h,ny-2));
                [nx,ny]=size(xm);
                xm=reshape(xm,1,nx*ny);
                ym=reshape(ym,1,nx*ny);

                % the constraint is the index of boundary elements

                X=[XX';[xm;ym]'];                    

                if str2double(vers(1:4))<2013
                    % olde code
                    fprintf('old')
                    dt=DelaunayTri(X,XC');
                    io=dt.inOutStatus;
                    p=dt.X';
                    t=dt.Triangulation';
                else
                     dt=delaunayTriangulation(X,XC');
                     p=dt.Points';  
                     t=dt.ConnectivityList';
                     io=dt.isInterior;                        
                end


                % Detect all triangles inside the object
                % Select all points from inside-triangles
                p=p(:,unique(t(:,io))');
                t=t(:,io);

                % New triangulation only with inside points
                % We need this iteration because of the  
                % numbering within the triangle-point relation


                if str2double(vers(1:4))<2013
                    % olde code                   
                    dt=DelaunayTri(p',XC');
                    io=dt.inOutStatus;
                    p=dt.X';
                    t=dt.Triangulation';
                else
                     dt=delaunayTriangulation(p',XC');
                     p=dt.Points';  
                     t=dt.ConnectivityList';
                     io=dt.isInterior;                        
                end
                obj.p=p;                    
                obj.e=XC;
                obj.e(3:4,:)=XE;
                % We have the same number of boundary segments as
                % arguments.
                obj.e(5,:)=NBDS;
                % Again, only interior triangles...
                obj.t=t(:,io);
                obj.t(4,:)=ones(1,size(obj.t,2)); 
                obj.jiggleMesh;
            end
        end
        
        function annulus(obj,r1,r2,h)  % HU
            fd=@(p) ddiff(dcircle(p,0,0,r2),dcircle(p,0,0,r1));
            pfix=[-1,-1;-1,1;1,-1;1,1];
%             fh=inline('min(sqrt( p(:,1).^2 + p(:,2).^2 ) , 1 )','p');
            fh=@(p) min(0.5+sqrt( p(:,1).^2 + p(:,2).^2 ) , 2 );
            [p,e,t]=distmesh2d(fd,fh,h,[-1,-1;1,1],pfix); 
             
            obj.p=p'; obj.t=t'; obj.e=e';
            
        end
        
        function holeInPlane(obj)
            % needs distmesh package
            %fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.4));
            fd=@(p) ddiff(dcircle(p,0,0,1),dcircle(p,0,0,0.4));
            pfix=[-1,-1;-1,1;1,-1;1,1];
%             fh=inline('min(sqrt( p(:,1).^2 + p(:,2).^2 ) , 1 )','p');
            fh=@(p) min(sqrt( p(:,1).^2 + p(:,2).^2 ) , 1 );
            [p,e,t]=distmesh2d(fd,fh,0.1,[-1,-1;1,1],pfix); 
             
            obj.p=p';
            obj.t=t';
            obj.e=e';
            
        end
        
        function circle(obj,r,h,x,y)
            % Creates the mesh for a circle geometry
            % obj.circle() Unitcircle with hmax=0.2
            % obj.circle(R) Circle with Radius=R and hmax=0.2
            % obj.circle(R,hmax)
            % obj.circle(R,hmax,xshist,yshift) 
            % Circle with center=(xshif,yshift)             
            % Code valid for MATLAB R>= 2013a by using new 
            % delaynayTriangulation class. For older Matlab
            % Releases it uses old DelaunayTri Class
            % (c) Uwe Prüfert            
              
            switch nargin
                case 1
                    h=0.2; 
                    R=linspace(0,1,round(1/h)+1);
                    x=0;
                    y=0;
                case 2   
                    h=r/5;                     
                    R=linspace(0,r,round(1/h)+1);
                    x=0;
                    y=0;
                case 3 
                    R=linspace(0,r,round(1/h)+1);
                    x=0;
                    y=0;
                case 5
                    R=linspace(0,r,round(1/h)+1);
                otherwise
                    obj.wrongNumberInputs.throw;
            end
            P=[];
            r=R(2)-R(1);
            for k=R 
                n=2*pi*k/r;
                s=linspace(0,2*pi,round(n));
                if length(s)>2
                    s(end)=[];
                end
                
                P=[P [k*sin(s);
                        k*cos(s)]];  
            end
            P =[P [0;0]];
            vers=version('-release');
            if str2double(vers(1:4)) >=2014
                dt=delaunayTriangulation(P');                        
                obj.p=dt.Points';
                obj.e=dt.freeBoundary';                
                obj.t=dt.ConnectivityList';                 
            else
                dt=DelaunayTri(P');                        
                obj.p=dt.X';
                obj.e=dt.freeBoundary';                
                obj.t=dt.Triangulation';  
                 
            end
            
            % Circle has only one boundary segment...
            % s runs from 0 to 1...
            
            
            s=linspace(0,1,obj.nEdges+1);
            obj.e(3,:)=s(1:end-1); 
            obj.e(4,:)=s(2:end);
            obj.e(5,:)=ones(1,size(obj.e,2));
            
            obj.t(4,:)=ones(1,size(obj.t,2));  
            obj.moveMesh(x,y);
        end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function square(obj,varargin)
% square(), square(h), square(xmin,xmax,ymin,ymax) square(xmin,xmax,ymin,ymax,h), 
% square(xmin,xmax,ymin,ymax,nx,ny)  % HU 
% No options => unitsqaure with h=0.1
% hmax given but no corner points => unitsquare with given hmax 
% left lower and right upper point given, but no hmax => square with hmax=0.1
% The boundary lenght parameter runs in anti clockwise order from 0 to 1 on everey boundary segment,
% independetly from its length!

% First handle all posible argument...     
gotnx=0; 
switch nargin
    case 1 % no options,unitsqaure with h=0.1
        xmin=0; xmax=1; ymin=0; ymax=1; h=0.1;                    
    case 2 % hmax but no coordinates, unitsquare with hmax                     
        xmin=0; xmax=1; ymin=0; ymax=1;  h=varargin{1};
    case 5 % left lower and right upper point given, but no hmax, square with hmax=0.1
        xmin=varargin{1}; xmax=varargin{2}; ymin=varargin{3}; ymax=varargin{4}; h=0.1; 
    case 6 % corners and hmax
        xmin=varargin{1}; xmax=varargin{2}; ymin=varargin{3}; ymax=varargin{4}; h=varargin{5};           
    case 7 % corners and nx,ny
        xmin=varargin{1}; xmax=varargin{2}; ymin=varargin{3}; ymax=varargin{4}; 
        nx=varargin{5}; ny=varargin{6}; gotnx=1; 
    otherwise
         obj.wrongNumberInputs.throw;
end
if ~gotnx; nx=max(2,round((xmax-xmin)/h)+1); ny=max(2,round((ymax-ymin)/h)+1); end         
xmesh=linspace(xmin,xmax,nx); ymesh=linspace(ymin,ymax,ny); [X,Y]=meshgrid(xmesh,ymesh);            
P=[reshape(X,1,nx*ny);reshape(Y,1,nx*ny)]; dt=delaunayTriangulation(P'); 
obj.p=dt.Points'; obj.e=dt.freeBoundary';  obj.t=dt.ConnectivityList'; 
% because we have an equidistant boundary, we can compute the boundary "s" in this simple way...
enx=linspace(0,1,nx); eny=linspace(0,1,ny); obj.e(3,1:nx-1)=enx(1:end-1); 
obj.e(3,nx:nx+ny-2)= eny(1:end-1); obj.e(3,nx+ny-1:nx+ny+nx-3)=fliplr(enx(1:end-1));
obj.e(3,nx+ny+nx-2:2*nx+2*ny-4)=fliplr(eny(1:end-1)) ;
obj.e(4,1:nx-1)=enx(2:end); obj.e(4,nx:nx+ny-2)= eny(2:end);
obj.e(4,nx+ny-1:nx+ny+nx-3)=fliplr(enx(2:end)); obj.e(4,nx+ny+nx-2:2*nx+2*ny-4)=fliplr(eny(2:end)) ;
% set the boundary segment number...
obj=setidssq(obj); 
obj.t(4,:)=ones(1,size(obj.t,2)); 
end 


function ccsquare(obj,varargin) % criss-cross-square, like square, but with rect-midpoints, HU
gotnx=0; 
switch nargin
    case 1 % no options,unitsqaure with h=0.1
        xmin=0; xmax=1; ymin=0; ymax=1; h=0.1;                    
    case 2 % hmax but no coordinates, unitsquare with hmax                     
        xmin=0; xmax=1; ymin=0; ymax=1;  h=varargin{1};
    case 5 % left lower and right upper point given, but no hmax, square with hmax=0.1
        xmin=varargin{1}; xmax=varargin{2}; ymin=varargin{3}; ymax=varargin{4}; h=0.1; 
    case 6 % corners and hmax
        xmin=varargin{1}; xmax=varargin{2}; ymin=varargin{3}; ymax=varargin{4}; h=varargin{5};           
    case 7 % corners and nx,ny
        xmin=varargin{1}; xmax=varargin{2}; ymin=varargin{3}; ymax=varargin{4}; 
        nx=varargin{5}; ny=varargin{6}; gotnx=1; 
    otherwise
         obj.wrongNumberInputs.throw;
end
if ~gotnx; nx=max(2,round((xmax-xmin)/h)+1); ny=max(2,round((ymax-ymin)/h)+1); end          
xmesh=linspace(xmin,xmax,nx); ymesh=linspace(ymin,ymax,ny); [X,Y]=meshgrid(xmesh,ymesh);            
P=[reshape(X,1,nx*ny);reshape(Y,1,nx*ny)];
dx=(xmax-xmin)/(nx-1); dy=(ymax-ymin)/(ny-1); % add rectangle midpoints to grid 
for i1=0:nx-2; for i2=0:ny-2; 
        xc=xmin+dx/2+i1*dx; yc=ymin+dy/2+i2*dy;
        cp=[xc;yc]; P=[P cp];
    end; end
P=P'; P=unique(P,'rows','stable'); P=P'; vers=version('-release');
if str2double(vers(1:4)) >=2014
    dt=delaunayTriangulation(P'); obj.p=dt.Points'; obj.e=dt.freeBoundary';  obj.t=dt.ConnectivityList';                 
else
    dt=DelaunayTri(P'); obj.p=dt.X'; obj.e=dt.freeBoundary'; obj.t=dt.Triangulation'; 
end
% because we have an equidistant boundary, we can compute the boundary "s" in this simple way...
enx=linspace(0,1,nx); eny=linspace(0,1,ny); obj.e(3,1:nx-1)=enx(1:end-1); 
obj.e(3,nx:nx+ny-2)= eny(1:end-1); obj.e(3,nx+ny-1:nx+ny+nx-3)=fliplr(enx(1:end-1));
obj.e(3,nx+ny+nx-2:2*nx+2*ny-4)=fliplr(eny(1:end-1)) ;
obj.e(4,1:nx-1)=enx(2:end); obj.e(4,nx:nx+ny-2)= eny(2:end);
obj.e(4,nx+ny-1:nx+ny+nx-3)=fliplr(enx(2:end)); obj.e(4,nx+ny+nx-2:2*nx+2*ny-4)=fliplr(eny(2:end)) ;

% set the boundary segment number...
obj.e(5,1:nx-1)=ones(1,nx-1); obj.e(5,nx:nx+ny-2)=2*ones(1,ny-1);
obj.e(5,nx+ny-1:nx+ny+nx-3)=3*ones(1,nx-1);
obj.e(5,nx+ny+nx-2:2*nx+2*ny-4)=4*ones(1,ny-1);           

obj.t(4,:)=ones(1,size(obj.t,2)); 
end     % ccsquare
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ellipse(obj)
            fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
            [p,e,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
            obj.p=p';
            obj.e=e';
            obj.t=t';
        end
        
        % turn and move actions
        function turnMesh(obj,arc)
            % turns mesh points 
            % obj.turnMesh(Arc)
            A=[cos(arc) -sin(arc);...
                sin(arc) cos(arc)];
            obj.p=A*obj.p;
        end
        
        function moveMesh(obj,x,y)
            % Moves meshpoints 
            % obj.moveMesh(xshift,yshift)
            if isscalar(x)&& isscalar(y)
                obj.p(1,:)=obj.p(1,:)+x;
                obj.p(2,:)=obj.p(2,:)+y;
            else
                obj.wrongFormat.throwAsCaller;
            end
        end
        
        % refine (and coarsening - later)
        
        function refineUniformly(obj,n)
            switch nargin
                case 1
                    n=1;
                case 2
                    if ~isscalar(n)
                        obj.wrongFormat.throw
                    end
                    
                otherwise
                    obj.wrongNumberInputs.throw
            end
            for k=1:n
                obj.refineMesh
            end
        end
        
        function refineMesh(obj,varargin)
            % Method for refinement of grid
            % obj.refineMesh();
            % obj.refineMesh(triangles_to_refine);            
 
            switch nargin
                case 1
                    toRefine=1:obj.nElements;
                case 2                    
                    toRefine=varargin{1};
                    if ~(isvector(toRefine)&&length(toRefine)<=obj.nElements)
                        % error if more elements are in refinement list as in
                        % the mesh
                        obj.wrongFormat.throw;
                    end
                otherwise
                    wrongNumberInputs.throw;
            end
             
            obj.refineRGB(toRefine);             
            obj.t=[obj.t; ones(1,size(obj.t,2))];            
        end
        
        function jiggleMesh(obj,prop,value)
            
            switch nargin
                case 1
                    obj.jiggle;
                case 3
                    switch prop
                        case 'n'
                            for k=1:value
                                obj.jiggle;
                            end
                        case 'quality'
                            k=0;
                            maxjiggle=15;
                            while min(getMeshQuality(obj))<value...
                                    &&k<maxjiggle
                               obj.jiggle;
                               k=k + 1;
                            end
                        otherwise
                            obj.wrongFormat.throwAsCaller;
                    end
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end
        
        % plot etc.
        
        function coneplot(obj,z0,a,varargin) % cone plot, HU                         
%ax=axescheck(varargin{:}); 
ax=newplot();  
x0=obj.p(1,:); y0=obj.p(2,:); r=sqrt(x0.^2+y0.^2); 
x=a*x0; y=a*y0; z=-r; 
patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
      'facevertexcdata',z0(:),'facecolor','interp', ...
        'edgecolor','none',  'parent',ax,varargin{2:end}); 
 set(gca,'Box','on'); colormap cool; view([-30 40]);   set(gca,'FontSize',14); 
 axis(ax,'tight'); axis(ax,'image'); 
end

function coneplotscal(obj,z0,a,s,varargin) % cone plot, HU                         
%ax=axescheck(varargin{:}); 
ax=newplot();  
x0=s*obj.p(1,:); y0=s*obj.p(2,:); r=sqrt(x0.^2+y0.^2); 
x=a*x0; y=a*y0; z=-r; 
patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
      'facevertexcdata',z0(:),'facecolor','interp', ...
        'edgecolor','none',  'parent',ax,varargin{2:end}); 
 set(gca,'Box','on'); colormap cool; view([-30 40]);   set(gca,'FontSize',14); 
 axis(ax,'tight'); axis(ax,'image'); 
end
        
        function plot(obj,varargin)
            % plot method for class gridd
            % obj.plot()
            % obj.plot(z)
            % obj.plot(arglist)
            % obj.plot(z,arglist)
            % arglist can contain pairs of parameters for controling plot,
            % see eg. help plot and the Matlab Reference page for plot.
            %               
       %     ax=axescheck(varargin{:});   
       ax=newplot(10);  
          
            if (nargin == 1) || (mod(nargin-1,2) == 0)               
                x=obj.p(1,:);
                y=obj.p(2,:);            
                z=zeros(size(x));
                patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
                    'facecolor',[1 1 1],...
                    'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                   'parent',ax,...
                    varargin{:}); 
                view(ax,2);
            else 
                x=obj.p(1,:); y=obj.p(2,:); z=varargin{1}; 
                if length(z) == obj.nPoints       
                    if 0
                    patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
                        'facevertexcdata',z(:),'facecolor','interp', ...
                        'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                        'parent',ax,varargin{2:end});
                    else % HU better scaling
                       
                     patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
                        'facevertexcdata',z(:),'facecolor','interp', ...
                        'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                        'parent',ax,varargin{2:end}); 
                    end
                    set(gca,'Box','on');               
                    view(3); 
                 %   grid(ax,'on');
                 %   colorbar;
                elseif length(z) == obj.nElements
                    x=reshape(obj.p(1,obj.t(1:3,:)),3,obj.nElements);
                    y=reshape(obj.p(2,obj.t(1:3,:)),3,obj.nElements);
                    z=[1;1;1]*z;
                    cz=z;
                    patch(x,y,z,cz,'lineStyle','none',varargin{2:end}); 
                    view(ax,3);
                    grid(ax,'on');
                    colorbar;
                end
            end  
           % axis(ax,'tight'); 
            %axis equal     
        end
        
        function toplot(obj,z0,R,rho,varargin) % parametric plot, HU                         
        ax=axescheck(varargin{:}); ax=newplot(ax);  
        x0=obj.p(1,:); y0=obj.p(2,:); x=(R+rho*cos(y0)).*cos(x0); 
        y=(R+rho*cos(y0)).*sin(x0);  z=rho*sin(y0); 
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','none', ... %get(ax,'DefaultSurfaceEdgeColor'),...
                'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; view([-30 40]);  set(gca,'FontSize',14); 
         axis(ax,'tight'); axis(ax,'image'); 
        end
        
        function spplot(obj,z0,R,varargin) % sphere plot, HU                         
        ax=axescheck(varargin{:}); ax=newplot(ax);  
        x0=obj.p(1,:); y0=obj.p(2,:);    
        x=R*cos(y0).*cos(x0); y=R*cos(y0).*sin(x0); z=R*sin(y0); 
        %size(x),size(y),size(z), size(z0)
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','none', ... %get(ax,'DefaultSurfaceEdgeColor'),...
                'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; view([-30 40]);   set(gca,'FontSize',14); 
         axis(ax,'tight'); axis(ax,'image'); 
        end
        function spplot0(obj,z0,R,varargin) % sphere plot, HU, black edges                     
        ax=axescheck(varargin{:}); ax=newplot(ax);  
        x0=obj.p(1,:); y0=obj.p(2,:);    
        x=R*cos(y0).*cos(x0); y=R*cos(y0).*sin(x0); z=R*sin(y0); 
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','k', ... %get(ax,'DefaultSurfaceEdgeColor'),...
                'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; view([-30 40]); set(gca,'FontSize',14); 
         axis(ax,'tight'); axis(ax,'image'); 
        end
        
        function spdplot(obj,z0,a,c,varargin) % ellipsoid plot, HU                         
        ax=axescheck(varargin{:}); ax=newplot(ax);  
        x0=obj.p(1,:); y0=obj.p(2,:); b=a; 
        x=a*cos(y0).*cos(x0); y=b*cos(y0).*sin(x0); z=c*sin(y0); 
        %size(x),size(y),size(z), size(z0)
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','none', ... %get(ax,'DefaultSurfaceEdgeColor'),...
                'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; view([-30 40]);   set(gca,'FontSize',14); 
         axis(ax,'tight'); axis(ax,'image'); 
        end
        
        function ellplot(obj,z0,a,b,c,varargin) % ellipsoid plot, HU                         
        ax=axescheck(varargin{:}); ax=newplot(ax);  
        x0=obj.p(1,:); y0=obj.p(2,:);    
        x=a*cos(y0).*cos(x0); y=b*cos(y0).*sin(x0); z=c*sin(y0); 
        %size(x),size(y),size(z), size(z0)
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','none', ... %get(ax,'DefaultSurfaceEdgeColor'),...
                'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; view([-30 40]);   set(gca,'FontSize',14); 
         axis(ax,'tight'); axis(ax,'image'); 
        end
        
        function cyplot(obj,z0,R,varargin) % cylinder plot, HU                         
        %ax=axescheck(varargin{:});     
        ax=newplot();  
        x0=obj.p(1,:); y0=obj.p(2,:);    
        x=R*cos(x0); y=R*sin(x0); z=y0; 
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','none'); %, ... %get(ax,'DefaultSurfaceEdgeColor'),...  'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; 
         view([30 40]);   set(gca,'FontSize',14); 
         if nargin>3; cmin=varargin{1};  cmax=varargin{2};  caxis([cmin cmax]); end 
         axis(ax,'tight'); axis(ax,'image'); 
        end
        
        function tocyplot(obj,z0,lz,R,varargin) % top of cylinder plot, HU                         
        %ax=axescheck(varargin{:});   
        ax=newplot();  
        x=R*obj.p(1,:); y=R*obj.p(2,:); z=lz*ones(1,length(y)); 
        patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
              'facevertexcdata',z0(:),'facecolor','interp', ...
                'edgecolor','none'); %, ... %get(ax,'DefaultSurfaceEdgeColor'),...                'parent',ax,varargin{2:end}); 
         set(gca,'Box','on'); colormap cool; 
         if nargin>3; cmin=varargin{1};  cmax=varargin{2};  caxis([cmin cmax]); end 
         view([30 40]);   set(gca,'FontSize',14); 
         axis(ax,'tight'); axis(ax,'image'); %colorbar
        end
        
       function [x,idx]=bdseg(obj,nSegment) % HU: mod of identifyBound. to extract Boundary segments 
        clf;  obj.plot; indx=find(obj.e(5,:)== nSegment);
        for k=1:length(obj.e(5,:))
            line(obj.p(1,obj.e(1:2,k)),obj.p(2,obj.e(1:2,k)),...
                'LineWidth',1,'color','blue')
        end 
        for k=1:length(indx)
            line(obj.p(1,obj.e(1:2,indx(k))),obj.p(2,obj.e(1:2,indx(k))),...
                'LineWidth',3,'color','red')
        end        
        idx=obj.e(1,indx(:)); x=obj.p(:,idx); 
        if nSegment<=max(obj.e(5,:))
            title(['Boundary segment no ',num2str(nSegment)]);
        else
            cla
            fprintf(['Sorry, there is not Boundary segment no ',num2str(nSegment),' in this geometry\n']);
        end
       end
     
        
        
        function identifyBoundarySegment(obj,nSegment)
            % obj.identifyBoundarySegment(nSegment)
            % plots the boundary segment number nSegment in red.
            % If the segment not exists, nothing will be highlighted.
            % (c) 2013 by Uwe Prüfert
            
            % Change log            
            % 2014/08/13 Fix a bug in identifyBoundarySegment when plotting
            % not ordered boundary segments
            % U.P.
            
            clf
            obj.plot;
            if nargin == 1 % on arguments
                for k=1:min(10,obj.nBoundarySegments)                    
                    figure()
                    obj.identifyBoundarySegment(k);   
                end
            else
                clf
                obj.plot;
                indx=find(obj.e(5,:) == nSegment);
    %             line(obj.p(1,obj.e(1:2,:)),obj.p(2,obj.e(1:2,:)));
                for k=1:length(obj.e(5,:))
                    line(obj.p(1,obj.e(1:2,k)),obj.p(2,obj.e(1:2,k)),...
                        'LineWidth',1,'color','blue')
                end 
                for k=1:length(indx)
                    line(obj.p(1,obj.e(1:2,indx(k))),obj.p(2,obj.e(1:2,indx(k))),...
                        'LineWidth',3,'color','red')
                end
                if nSegment<=max(obj.e(5,:))
                    title(['Boundary segment no ',num2str(nSegment)]);
                else
                    cla
                    fprintf(['Sorry, there is not Boundary segment no ',num2str(nSegment),' in this geometry\n']);
                end
            end
        end
        
        
        function quiver(obj,dx,dy)
%             nt=size(obj.t,2);
                    x=reshape(obj.p(1,obj.t([1 2 1 3 2 3],:)),6,obj.nElements);              
                    z=zeros(size(x));
                    y=reshape(obj.p(2,obj.t([1 2 1 3 2 3],:)),6,obj.nElements);  
                    h=patch(x,y,z);
                    set(h,'FaceColor','none','EdgeColor',[0.9 0.9 0.9]);
                    view(0,90)
            for k=1:obj.nPoints
                line([obj.p(1,k),obj.p(1,k)+0.05*dx(k)],...
                    [obj.p(2,k),obj.p(2,k)+0.05*dy(k)],'Color','r');
            end
        end
        
        % utilities
        
        function boundaryPoints=getBoundaryPoints(obj)
            %boundaryPoints=obj.getBoundaryPoints
            % Gives back the point on the boundary of a grid.
            % (c) 2013 by Uwe Prüfert
            boundaryPoints=obj.p(:,unique(obj.e(1:2,:)));        
        end
        
        function n=nBoundaryPoints(obj)
            n=length(obj.getBoundaryPoints);
        end
        
        function index= getBoundaryPointsIndex(obj)
            %index= obj.getBoundaryPointsIndex
            % Gives back the index of the points on the boundary of a grid.
            % (c) 2013 by Uwe Prüfert
            index=unique(obj.e(1:2,:))';        
        end
        
        function index= getInnerPointsIndex(obj)
            %index=obj.getInnerPointsIndex
            % Gives back the index of the points not on the boundary of a grid.
            % (c) 2013 by Uwe Prüfert
             
            index=1:obj.nPoints;
            index(obj.getBoundaryPointsIndex)=[];      
        end
        
        function innerPoints=getInnerPoints(obj)
            %innerPoints=obj.getInnerPoints
            % Gives back the  points not the boundary of a grid.
            % (c) 2013 by Uwe Prüfert
            indx=obj.getInnerPointsIndex;
            if ~isempty(indx)
                innerPoints=obj.p(:,indx);
            else
                warning('grid2D:emptyInnerPointsIndex',...
                    'Grid contains no inner points.');
            end
        end
        
        function val=getMeshQuality(obj)
            % angle=getInnerAngles(obj)
            
            % a,b,c are the square of the edge lenghts
            a=sum((obj.p(:,obj.t(2,:))-obj.p(:,obj.t(1,:))).^2);
            b=sum((obj.p(:,obj.t(3,:))-obj.p(:,obj.t(2,:))).^2);
            c=sum((obj.p(:,obj.t(3,:))-obj.p(:,obj.t(1,:))).^2);
           
            val=sum([(b+c-a)./(2*sqrt(b.*c));...
                (a+c-b)./(2*sqrt(a.*c));...
                (a+b-c)./(2*sqrt(a.*b))]);
            
        end
        
        % coefficients 
        
        function [bvalvec]=convCoefficientsMpt(obj,b)
            %computes the value of the coefficient b in the center of every triangle
            % coefficent can be a vector of dim 2 x 1, or a cell array of dim 2 x 1 
            
            % 'symbolic' variables x, y  and t are neccesary for evaluation of string 
            % objects like c='sin(x)' etc.
            p=obj.p;
            t=obj.t;
            
            x=p(1,:);
            y=p(2,:);
            
            % number of points and triangles
%             n=length(p(2,:));
%             nt=length(t(1,:));
            % check the class and  size of b
            
            if ~((max(size(b))>=2)&&(min(size(b))>=1))
                ME=MException('ccoefficients:wrongCoefficientDefinition',...
                    ' b must be a vector');
                throw(ME);
            end
            % two cases:
            % cell-array - entries are strings, or doubles
            if isa(b,'cell')
                b1=b{1};
                b2=b{2};
                
            elseif isa(b,'double')
                b1=b(1,:);
                b2=b(2,:);
            else
                ME=MException('ccoefficients:wrongCoefficientDefinition',...
                    ' Wrong coefficient definition');
                throw(ME);
            end
            
            if isa(b1,'function_handle') || isa(b1,'inline')
                bval=feval(b1,x,y);
            elseif isa(b1,'char'),
                bval=eval(b1).*ones(1,obj.nPoints);
            elseif isa(b1,'numeric')
                if length(b1)==obj.nPoints
                    % c vektor and defined in p
                    bval=b1;
                elseif length(b1)==1,
                    % skalar
                    bval=b1*ones(obj.nElements,1);
                elseif length(b1)==obj.nElements,
                    bval=b1;
                else
                    MException('ccoefficients:wrongSize',...
                        'wrong sized b(1)').throwAsCaller;
                     
                end
            elseif isa(b1,'inline')
                bval=b(x,y);
            else
                MException('ccoefficients:wrongSize',...
                    'wrong formated b(1)').throwAsCaller;
                 
            end
            if length(bval) == obj.nElements
                % b(1) is a vektor and defined in center of mass of triagle
            else
                dimb=size(bval);
                if dimb(1) == 1
                    bval=bval';
                end
                bval=obj.point2Center(bval);
            end
            dimb=size(bval);
            if dimb(1) == 1
                bval=bval';
            end
            % first column of bvalvec
            bvalvec=bval;
            
            if isa(b2,'function_handle') || isa(b2,'inline')
                bval=feval(b2,x,y);
            elseif isa(b2,'char'),
                bval=eval(b2).*ones(1,obj.nElements);
            elseif isa(b2,'numeric')
                if length(b2)==obj.nPoints
                    % c vektor and defined in p
                    bval=b2;
                elseif length(b2)==1,
                    % skalar
                    bval=b2*ones(obj.nElements,1);
                elseif length(b2)==obj.nElements
                    bval=b2;
                else
                    ME=MException('ccoefficients:wrongSize','wrong sized b(1)');
                    throw(ME);
                end
            elseif isa(c,'inline')
                bval=b(x,y);
            else
                ME=MException('ccoefficients:wrongSize','wrong formated b(1)');
                throw(ME);
            end
            
            if length(bval) == obj.nElements
                % b(1) is a vektor and defined in center of mass of triagle
                
            else
                dimb=size(bval);
                if dimb(1) == 1
                    bval=bval';
                end
                bval=obj.point2Center(bval);
            end
            dimb=size(bval);
            if dimb(1) == 1
                bval=bval';
            end
            bvalvec=[bvalvec,bval];
        end
        
        function [cval,aval,fval]=aCoefficientsMpt(obj,c,a,f)
            %computes the value of the coefficients in the center of every triangle
            % 'symbolic' variables x, y  are neccesary for evaluation of string 
            % objects like c='sin(x)' etc.           
            p=obj.p;    t=obj.t;
            switch class(c)
                case 'function_handle'
                    midp=obj.midpts;  x=midp(1,:); y=midp(2,:); cval=feval(c,x,y);
                    [rows,cols]=size(cval);
                    if rows == 1;   cval=[cval;cval]; end % scalar
                    if max(rows,cols)==obj.nPoints
                        cval=[obj.point2Center(cval(1,:)); obj.point2Center(cval(2,:))];
                    end                      
              case 'double' % four cases: (i) scalar (ii) 2 x 2 matrix 
                             %  (iii) vector length np, (iv) vector length nt
                    [rows,cols]=size(c); 
                    switch max(rows,cols)
                        case 1; cval=c(ones(2,obj.nElements));
                        case 2; % constant matrix 
                    if c(1,2)==0; % diagonal case 
                     c1=c(1,1); c2=c(2,2); %pause
                     cval(1,:)=c1(ones(1,obj.nElements));                     
                     cval(2,:)=c2(ones(1,obj.nElements)); 
                    else % 
                     c1=c(1,1); c2=c(1,2); c3=c(2,1); c4=c(2,2); 
                     cval(1,:)=c1(ones(1,obj.nElements));  % pa_x(c1*u_x)                 
                     cval(2,:)=c2(ones(1,obj.nElements));  % pa_x(c2*u_y)
                     cval(3,:)=c3(ones(1,obj.nElements)); % pa_y(c3*u_x)
                     cval(4,:)=c4(ones(1,obj.nElements));  % pa_y(c4*u_y)                            
                    end 
                        case obj.nPoints
                            cval=obj.point2Center(c); cval=[cval;cval];
                        case obj.nElements  % the good one, nothing to do
                            cval=[c(:)';c(:)'];
                          otherwise % HU: matrix, [c1(1:np c2(1:np); c3(1:np) c4(1:np)] 
                    %               or  [c1(1:nt c2(1:nt); c3(1:nt) c4(1:nt)] 
                    np=obj.nPoints; nt=obj.nElements; 
                    if cols==2*np 
                        cval(1,:)=obj.point2Center(c(1,1:np));
                        cval(2,:)=obj.point2Center(c(1,np+1:2*np));
                        cval(3,:)=obj.point2Center(c(2,1:np));
                        cval(4,:)=obj.point2Center(c(2,np+1:2*np));
                   elseif rows==2*np; % diagonal matrix
                        cval(1,:)=obj.point2Center(c(1:np));
                        cval(2,:)=zeros(1,nt); 
                        cval(3,:)=zeros(1,nt); 
                        cval(4,:)=obj.point2Center(c(np+1:2*np));
                    else 
                        cval(1,:)=c(1,1:nt); cval(2,:)=c(1,nt+1:2*nt);
                        cval(3,:)=c(2,1:nt); cval(4,:)=c(2,nt+1:2*nt);
                    end
                    end
                case 'char'
                    % must be a single char symbolizing
                    % the coefficent function.
                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".  
                    % Hence, the next warning can be ignored!
                    midp=obj.midpts;
                    x=midp(1,:);
                    y=midp(2,:); %#ok<*NASGU>
                    try
                        cval=eval(c);
                    catch ME
                        throw(ME);
                    end
                    [rows,cols]=size(cval);
                    switch max(rows,cols)
                        case 1
                            % c is a constant like 'pi'
                            cval=cval(ones(2,obj.nElements));
                        case obj.nElements
                            cval=[cval;cval];
                        otherwise
                            throw(obj.wrongFormat);
                    end
                case 'cell'
                    error('Sorry, Cell array input not jet implemented!')
                otherwise
                    throw(obj.wrongClass)
            end
            % repeat code from case 'c'...
            switch class(a)
                case 'function_handle'
                    midp=obj.midpts;
                    x=midp(1,:);
                    y=midp(2,:);
                    aval=feval(a,x,y);
                    [rows,cols]=size(aval);
                    if max(rows,cols)==obj.nPoints
                        aval= obj.point2Center(aval);
                        aval=aval(:)';
                    end  
                case 'double'
                    [rows,cols]=size(a);
                    switch max(rows,cols)
                        case 1
                            aval=a(ones(1,obj.nElements));            
                        case obj.nPoints
                            aval=obj.point2Center(a);
                            aval=aval(:)';
                            %                            
                        case obj.nElements
                            % the good one, nothing to do
                            aval=a(:)';
                        otherwise
                            obj.wrongFormat.throwAsCaller
                    end
                case 'char'         
                    midp=obj.midpts;
                    x=midp(1,:);
                    y=midp(2,:);
                    try
                        aval=eval(a);
                    catch ME
                         ME.throwAsCaller;
                    end
                    [rows,cols]=size(aval);
                    switch max(rows,cols)
                        case 1
                            aval=aval(ones(1,obj.nElements));             
                        case obj.nElements
                            % work already done 
                        otherwise
                            finiteElements.wrongFormat.throwAsCaller;
                    end
                case 'cell'
                    obj.wrongClass.throwAsCaller;
                otherwise
                    obj.wrongClass.throwAsCaller;
            end 
           
            switch class(f)
                case 'function_handle'
                    midp=obj.midpts;
                    x=midp(1,:);
                    y=midp(2,:);
                    fval=feval(f,x,y);
                    [rows,cols]=size(fval);
                    if max(rows,cols)==obj.nPoints
                        fval= obj.point2Center(fval);
                    end  
                case 'double'  
                    [rows,cols]=size(f);          
                    switch max(rows,cols)
                        case 1
                            fval=f(ones(1,obj.nElements));
                                                          
                        case obj.nPoints
                            fval=obj.point2Center(f);  
                                                   
                        case obj.nElements
                            % the good one, nothing to do
                            fval=f;         
                        otherwise
                            obj.wrongFormat.throw
                    end
%                     max(fval)
                case 'char'        
                    midp=obj.midpts;
                    x=midp(1,:);
                    y=midp(2,:);
                    try
                        fval=eval(f);
                    catch ME
                        throw(ME);
                    end

                    [rows,cols]=size(fval);
                    switch max(rows,cols)
                        case 1
                           fval=fval(ones(1,obj.nElements));
                        case obj.nElements
                            % work already done
                        otherwise
                            throw(obj.wrongFormat);
                    end
                case 'cell'
                    error('Sorry, Cell array input not jet implemented!')
                otherwise
                    throw(obj.wrongClass)
            end

        end    
        
        % utilities: diameter, center of triangle, has triangle a boundary edge etc.
        
        function diam=triangleDiameters(obj)
            % diam=obj.triangleDiameters() 
            % Computes the (outer) diameter of all triangles
            d1=(obj.p(:,obj.t(1,:))-obj.p(:,obj.t(2,:)));
            d2=(obj.p(:,obj.t(1,:))-obj.p(:,obj.t(3,:)));
            d3=(obj.p(:,obj.t(2,:))-obj.p(:,obj.t(3,:)));
            diam(1,:)=sqrt(diag(d1'*d1));
            diam(2,:)=sqrt(diag(d2'*d2));      
            diam(3,:)=sqrt(diag(d3'*d3));             
            diam=max(diam);
        end 
        
        function mpt=midpts(obj)
            % mpt=obj.midpts()
            % helper method of class grid2D
            % computes the coordinates of the midpoints
            % of triangles of gt             
            mpt=1/3*(obj.p(:,obj.t(1,:))...
                +obj.p(:,obj.t(2,:))...
                +obj.p(:,obj.t(3,:)));
        end            
        
        function b=isBoundaryTriangle(obj)
            %  b=isBoundaryTriangle(gt) 
            % Method that indicates the boundary points.
            % b is a vector of length(#triangles)    %
                        
            % Algorithmus: falls in der k-ten Spalte von N 
            % Nullen stehen, hat das Dreieck weniger als drei
            % Nachbarn, ist also ein Dreieck mit Randkante            
            N=obj.neighbours();            
            b=false(1,size(obj.t,2));
            for k=1:size(N,2)
                n=0; % keine Randkanten
                for j=1:3
                    if N(j,k)==0,
                        n=n+1;
                    end
                end
                b(k)=(n>0);
            end
        end    
        
        function indx=pointToElementIndex(obj,pt)
            % pointToElementIndex gives back the index of the triangle
            % containing the point pt pt must be a vector of length
            % two.
            % indx=isPointInTriangle(obj,pt)
           
            b= obj.isPointInTriangle(pt);
            indx=find(b);            
        end
        
        function b=isPointInDomain(obj,pt)
            % isPointInDomain - true if the point pt is in at least
            % one triangle of the triangulation
            % 
            b=any(obj.isPointInTriangle(pt));
        end
        
        
        
        
        function varargout= gradient(obj,f,varargin)  
            % gradient -- computes the gradient of f
            % obj.gradient(f)
            % f must be a vector of length number of points
            
            % Default values: 100 Points for meshgrid and evaluation on the
            % grids points.
            n=100;
            p1=obj.p(1,:);
            p2=obj.p(2,:);
            
            % insert here options
            
            if nargout == 0
                % without return argument we can skip everything ;-)
                % maybe in future throw here an error
                return
            end
            
            f=full(f);
            
            vers=version('-release');
            if any(vers(1:4) >= '2014')                 
                % compute interpolant for data "f"
                F=scatteredInterpolant(obj.p(1,:)',obj.p(2,:)',f);
                % meshgrid
                tx=linspace(min(obj.p(1,:)),max(obj.p(1,:)),n)';
                ty= linspace(min(obj.p(2,:)),max(obj.p(2,:)),n)';
                [X,Y]=meshgrid(tx,ty); 
                % Evaluate F on new grid
                z=F(X,Y);
                [dx,dy]=gradient(z,tx,ty);
                
                % COmpute interpolant for the gradeint
                DX=scatteredInterpolant(reshape(X,n^2,1),...
                                        reshape(Y,n^2,1),...
                                        reshape(dx,n^2,1));

                DY=scatteredInterpolant(reshape(X,n^2,1),...
                                        reshape(Y,n^2,1),...
                                        reshape(dy,n^2,1));
                % Evaluate interpoalnts
                dx=DX(p1,p2);
                dy=DY(p1,p2);
            else            
                F=TriScatteredInterp(obj.p(1,:)',obj.p(2,:)',f);
                tx=linspace(min(obj.p(1,:)),max(obj.p(1,:)),n)';
                ty= linspace(min(obj.p(2,:)),max(obj.p(2,:)),n)';

                         
                [X,Y]=meshgrid(tx,ty);            
                z=F(X,Y);           
                [dx,dy]=gradient(z,tx,ty);   

                DX=TriScatteredInterp(reshape(X,n^2,1),...
                                        reshape(Y,n^2,1),...
                                        reshape(dx,n^2,1));

                DY=TriScatteredInterp(reshape(X,n^2,1),...
                                        reshape(Y,n^2,1),...
                                        reshape(dy,n^2,1));

                dx=DX(obj.p(1,:),obj.p(2,:));
                dy=DY(obj.p(1,:),obj.p(2,:));
            end
           
            switch nargout
                case 0
                    % do nothing                    
                case 1                   
                    varargout{1}=[dx;dy];
                case 2
                    varargout{1}=dx;
                    varargout{2}=dy;
                otherwise
            end
        end  
        function [sidelength,area]=sideLengthAndArea(obj)
        %Method to compute  side lengths and areas of triangles.
            dx=zeros(3,obj.nElements);
            dy=zeros(3,obj.nElements);
            sidelength=zeros(3,obj.nElements);
            for k=1:3,
                k1=rem(k,3)+1;
                k2=rem(k1,3)+1;
                dx(k,:)=obj.p(1,obj.t(k1,:)) - obj.p(1,obj.t(k2,:));
                dy(k,:)=obj.p(2,obj.t(k1,:)) - obj.p(2,obj.t(k2,:));
                sidelength(k,:)=sqrt(dx(k,:).*dx(k,:) + dy(k,:).*dy(k,:));
            end;
            area=0.5*abs(dx(1,:).*dy(2,:) - dx(2,:).*dy(1,:));
        end
    end
    
    methods(Access=protected)
        function neighbourIndex=getNeighbourPointsIndex(obj,k,d)
            %neighbourIndex=obj.getNeighbourPointsIndex(k)
            % Computes the index of all neighboured points of  
            % point k.            
            %neighbourIndex=obj.getNeighbourPointsIndex(k,d)
            % Computes the index of all neighboured points of  
            % point k in Domain d.
            
            [~,indx]=find(obj.t(1:3,:)==k);
            if (nargin == 3)
                indx=indx(obj.t(4,indx) == d);
            end
            neighbourIndex=unique(obj.t(1:3,indx))';
            neighbourIndex(neighbourIndex==k)=[];
        end
    end 
    
    methods(Access=private)
        
        function jiggle(obj)      
            %obj.jiggle
            % Jiggles the mesh
            % (c) 2013 by Uwe Prüfert
             
            innerPoints=obj.getInnerPointsIndex;
            for k=innerPoints                
                indx=obj.getNeighbourPointsIndex(k);
                % mean value 
                obj.p(:,k)=sum(obj.p(:,indx),2)/length(indx);
            end
            
        end
        
        
        
        
        function b=isPointInTriangle(obj,pt)
            % isPointInTriangle -  gives back a bool vector of length nTri.
            % b ist true if the point pt is i the triangle of a given
            % mesh object.
            % b=obj.isPointInTriangle(pt) 
            % Algortihm: Transform the point wrt the transformation of 
            % the triangle into the unit triangle  and
            % decide based on the relation pt_x < 0 pt_y < 0 ... etc.
            p=obj.p;  t=obj.t;   p1=p(:,(t(1,:)));
            p2=p(:,(t(2,:)));
            p3=p(:,(t(3,:)));
            x21=p2(1,:)-p1(1,:);
            x31=p3(1,:)-p1(1,:);
            y21=p2(2,:)-p1(2,:);
            y31=p3(2,:)-p1(2,:);            
            J=x21.*y31-x31.*y21; 
            % transform the point into the unit-triangle and
            % decide
            ptx=pt(1)-p1(1,:);
            pty=pt(2)-p1(2,:);
            ptux=1./J.*(y31.*ptx-x31.*pty);
            ptuy=1./J.*(x21.*pty-y21.*ptx);
            b=~(sum([ptux;ptuy])>1 | ptux<0 | ptuy<0);
        end
        
        function N=neighbours(obj,varargin)
            % N=neigbours(gt[,indx]) 
            % computes to every triangle the 
            % index of neighbored triangles              
%             nt=length(obj.t(1,:));            
            if nargin==2,
                indx=varargin{1};
            else
                indx=1:obj.nElements;
            end            
            N=sparse_null(3,obj.nElements);
            for k=1:length(indx),
                nk=obj.ent(indx(k));
                [i]=find(nk==k);
                nk(i)=[];
                N(1:length(nk),indx(k))=nk;
            end
        end
        
        function intl=ent(obj,it)
            % ent Indices of triangles neighboring
            % a given set of triangles. 
%             nt=size(obj.t,2);
            obj.nPoints=size(obj.p,2);
            switch nargin
                case 1
                    % all neighbours
                    it=1:obj.nElements;
                case 2
                    % okay
                otherwise
                    obj.wrongNumberInputs.throw;
            end   
            it1=ones(1,obj.nElements);
            it1(it)=zeros(size(it));
            it1=find(it1);  
            ip1= obj.t(1,it)';
            ip2= obj.t(2,it)';
            ip3= obj.t(3,it)';
            
            % Make a connectivity matrix.
            A=sparse(ip1,ip2,1,obj.nPoints,obj.nPoints);
            A=A+sparse(ip2,ip3,1,obj.nPoints,obj.nPoints);
            A=A+sparse(ip3,ip1,1,obj.nPoints,obj.nPoints);
            
            ntl=zeros(1,nnz(A-A')); % a slight overestimate
            nnt=0;          
            
            for i=1:length(it1),
                if A( obj.t(2,it1(i)), obj.t(1,it1(i))) || ...
                        A( obj.t(3,it1(i)), obj.t(2,it1(i))) || ...
                        A( obj.t(1,it1(i)), obj.t(3,it1(i)))
                    nnt=nnt+1;
                    ntl(nnt)=it1(i);
                end
            end
            intl=sort([it ntl(1:nnt)]);
        end 
        
        function  refineRGB(obj,markedElements)
            %refineRGB: local refinement of finite element mesh by red-green-blue
            %           refinement, where marked elements are red-refined.
            %
            %Usage:
            %
            %obj.refineRGB(markedElements)
            %  
            %    This methods bases on a code taken from the paper 
            %    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
            %    by S. Funken, D. Praetorius, and P. Wissgott. 
            % Original Authors: 
            % 
            %    S. Funken, D. Praetorius, P. Wissgott  10-07-08
            %  Adaption for grid2D class by Uwe Pruefert   2013
         
            coordinates=obj.p';
            boundary=obj.e';
            elements=obj.t(1:3,:)';
             
            nE=size(elements,1);
            %*** Sort elements such that first edge is longest
            dx=coordinates(elements(:,[2,3,1]),1)-coordinates(elements,1);
            dy=coordinates(elements(:,[2,3,1]),2)-coordinates(elements,2);
            [~,idxMax]=max(reshape(dx.^2+dy.^2,nE,3),[],2);
            idx=( idxMax==2 );
            elements(idx,:)=elements(idx,[2,3,1]); 
            idx=( idxMax==3 );
            elements(idx,:)=elements(idx,[3,1,2]); 
            %*** Obtain geometric information on edges
             
            [edge2nodes,element2edges,boundary2edges{1}] ...
               =provideGeometricData(elements,boundary(:,1:2));
             
            %*** Mark edges for refinement
            edge2newNode=zeros(max(max(element2edges)),1);
            try; rlong=obj.rlong; catch; rlong=0; end 
            if rlong;  edge2newNode(element2edges(markedElements,1))=1;
            else edge2newNode(element2edges(markedElements,:))=1;
            end        
            swap=1; %swap=[]; pause % HU
            while ~isempty(swap)
                markedEdge=edge2newNode(element2edges);
                swap=find( ~markedEdge(:,1) & (markedEdge(:,2) | markedEdge(:,3)) );                
                edge2newNode(element2edges(swap,1))=1;
            end
            %*** Generate new nodes
        
            edge2newNode((edge2newNode)>0)=size(coordinates,1) + ...
                (1:nnz(edge2newNode));
            idx=find(edge2newNode);
            coordinates(edge2newNode(idx),:) ...
               =(coordinates(edge2nodes(idx,1),:)+...
                coordinates(edge2nodes(idx,2),:))/2;
            %*** Refine boundary elements 
            % In constrast to the original code,
            % we must extend the boundary matrix
            % by adding informations on the number
            % of boundary and the start and end
            % parameter. Importand for evaluating
            % u=f(s) statements... U.P.
            if ~isempty(boundary)
                newNodes=edge2newNode(boundary2edges{1});  
                markedEdges=find(newNodes);
                % We compute the   midpoint of
                % edge. It will be the end point
                % of the first end the starting
                % point of the second set of new
                % edges. U.P.
                newEdgeMidpoint=(boundary(markedEdges,4)+boundary(markedEdges,3))/2;
                
                if ~isempty(markedEdges)
                    boundary=[boundary(~newNodes,:); ...
                                boundary(markedEdges,1),newNodes(markedEdges),boundary(markedEdges,3),newEdgeMidpoint,boundary(markedEdges,5:end); ...
                                newNodes(markedEdges),boundary(markedEdges,2),newEdgeMidpoint,boundary(markedEdges,4:end)]; 
                end
            end
            newNodes=edge2newNode(element2edges);
            %*** Determine type of refinement for each element
            markedEdges=(newNodes~=0);
            none=~markedEdges(:,1);
            bisec1 =( markedEdges(:,1) & ~markedEdges(:,2) & ~markedEdges(:,3) );
            bisec12=( markedEdges(:,1) &  markedEdges(:,2) & ~markedEdges(:,3) );
            bisec13=( markedEdges(:,1) & ~markedEdges(:,2) &  markedEdges(:,3) );
            red    =( markedEdges(:,1) &  markedEdges(:,2) &  markedEdges(:,3) );
            %*** Generate element numbering for refined mesh
            idx=ones(nE,1);
            idx(bisec1) =2; %*** green=newest vertex bisection of 1st edge
            idx(bisec12)=3; %*** blue (right)=newest vertex bisection of 1st and 2nd edge
            idx(bisec13)=3; %*** blue (left)=newest vertex bisection of 1st and 3rd edge
            idx(red)    =4; %*** red refinement
            idx=[1;1+cumsum(idx)];
            %*** Generate new elements
            newElements=zeros(idx(end)-1,3);
            newElements(idx(none),:)=elements(none,:);
            newElements([idx(bisec1),1+idx(bisec1)],:) ...
               =[elements(bisec1,3),elements(bisec1,1),newNodes(bisec1,1); ...
                   elements(bisec1,2),elements(bisec1,3),newNodes(bisec1,1)];
            newElements([idx(bisec12),1+idx(bisec12),2+idx(bisec12)],:) ...
               =[elements(bisec12,3),elements(bisec12,1),newNodes(bisec12,1); ...
                   newNodes(bisec12,1),elements(bisec12,2),newNodes(bisec12,2); ...
                   elements(bisec12,3),newNodes(bisec12,1),newNodes(bisec12,2)]; 
            newElements([idx(bisec13),1+idx(bisec13),2+idx(bisec13)],:) ...
               =[newNodes(bisec13,1),elements(bisec13,3),newNodes(bisec13,3); ...
                   elements(bisec13,1),newNodes(bisec13,1),newNodes(bisec13,3); ...
                   elements(bisec13,2),elements(bisec13,3),newNodes(bisec13,1)];
            newElements([idx(red),1+idx(red),2+idx(red),3+idx(red)],:) ...
               =[elements(red,1),newNodes(red,1),newNodes(red,3); ...
                   newNodes(red,1),elements(red,2),newNodes(red,2); ...
                   newNodes(red,3),newNodes(red,2),elements(red,3); ...
                   newNodes(red,2),newNodes(red,3),newNodes(red,1)];

            obj.p=coordinates';
            obj.e=boundary';
            obj.t=newElements'; 
            
            % local function
            function [edge2nodes,element2edges,boundaries] ...
                   =provideGeometricData(elements,boundaryG)
                %provideGeometricData: returns geometric data for finite element mesh
                %
                %Usage:
                %
                %[edges2nodes,element2edges,dirichlet2edges,neumann2edges] ...
                %   =provideGeometricData(elements,dirichlet,neumann)
                %
                %Authors:
                % 
                %    S. Funken, D. Praetorius, P. Wissgott  10-07-08
          
                nEE=size(elements,1);
                nB=nargin-1;
                %*** Node vectors of all edges (interior edges appear twice) 
                I=elements(:); 
                J=reshape(elements(:,[2,3,1]),3*nEE,1); 
                %*** Symmetrize I and J (so far boundary edges appear only once)
                pointer=[1,3*nEE,zeros(1,nB)];
                 
                 
                I=[I;boundaryG(:,2)];  
                J=[J;boundaryG(:,1)]; 
                     
                pointer(3)=pointer(2) + size(boundaryG,1);
               
                %*** Create numbering of edges
                idxIJ=find(I < J); 
                edgeNumber=zeros(length(I),1);
                edgeNumber(idxIJ)=1:length(idxIJ); 
                idxJI=find(I > J);
                number2edges=sparse(I(idxIJ),J(idxIJ),1:length(idxIJ));
                [foo{1:2},numberingIJ]=find(number2edges);
                [foo{1:2},idxJI2IJ]=find(sparse(J(idxJI),I(idxJI),idxJI) );  
               % idxJI2IJ, numberingIJ
                edgeNumber(idxJI2IJ)=numberingIJ;
                %*** Provide element2edges and edge2nodes
                element2edges=reshape(edgeNumber(1:3*nEE),nEE,3);
                edge2nodes=[I(idxIJ),J(idxIJ)];
                %*** Provide boundary2edges
                
                boundaries=edgeNumber(pointer(2)+1:pointer(3));
                
            end
        end
    end
    
    methods(Static)  
        % here bcoefficient is STATIC         
        function[qval,gval,hval,rval]=boundCoefficients(p,b)
            % compute the boundary coefficients
            % 
            % q,g,h,r are everything that can be evaluated by eval of feval, The
            % independent variables must be named by x,y (Euklidian) or by s (arc-length
            % parametrization)

           if length(p) == 2,           
                x=p(1);
                y=p(2);
            else
                s=p;
            end
            m=b(2);
            qval=0;
            gval=0;
            hval=0;
            rval=0;
            lengthq=b(3);
            lengthg=b(4);
            if m == 0 % only Neumann BCs
                qval=eval(char(b(5:5+lengthq-1)));
                gval=eval(char(b(5+lengthq:5+lengthq+lengthg-1)));
            else % only Dirichlet BCs
                lengthh=b(5);
                lengthr=b(6);
                char(b(9:9+lengthh-1));
                hval=eval(char(b(9:9+lengthh-1)));
                rval=eval(char(b(9+lengthh:9+lengthh+lengthr-1)));
            end  
        end         
    end

    % we don't need this stuff, so we better hide it
    methods(Hidden)
        % all methods from handle
        % intended not to shown...
       %  addlistener(obj)
       %  ge(obj)
       %  le(obj)
       %  ne(obj)
       %  notify(obj)
       %  gt(obj)
       %  eq(obj)
       %  findobj(obj)
       %  findprop(obj)
       %  lt(obj)        
    end  
end
