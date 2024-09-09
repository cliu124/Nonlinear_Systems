classdef grid3D < gridd & plotUtils3D
%grid3D 3D grid for used with FEM-OO Toolbox
% Class for tetrahedral grids.   
% Change log:
% 01/2014 Added methods: cutawyPlot, plotIso, and plotSlices
% 02/2014 Added methods:  refineMesh, redRefine and provideEdgeData.
%         plotFaces: Now interpolated insted of flat data
%         plotting.
%         Added methods: ficheraCube
%

methods(Access=public)       
    function obj=grid3D()
        % constructor method
        % only empty call allowed
        % obj=grid3D()

        % only  in geometry
        obj.nPointsInElements=4;
        switch nargin
            case 0
                % empty object
            otherwise
                obj.wrongNumberInputs.throw;
        end
    end 
        
function pts2dom(obj,P) % point-cloud to PTE 
dt=delaunayTriangulation(P'); obj.p=dt.Points'; obj.e=dt.freeBoundary'; 
obj.t=dt.ConnectivityList'; 
obj.t(5,:)=ones(1,size(obj.t,2));  obj=setids1(obj); 
end         
        
function bar(obj,varargin)           
% method that meshes a 3D bar.
% obj.bar(xmin,xmax,ymin,ymax,zmin,zmax)
% obj.bar(xmin,xmax,ymin,ymax,zmin,zmax,hmax)
% obj.bar(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz)  % HU 
% (c) 2013 by Uwe Prüfert
gotnx=0; 
switch nargin
    case {7,8,10}
        switch nargin
            case 7; hmax=0.1;
            case 8; hmax=varargin{7}; 
            case 10; nx=varargin{7}; ny=varargin{8}; nz=varargin{9};gotnx=1; 
        end
        if ~gotnx
        nx=max(2,ceil((varargin{2}-varargin{1})/hmax)+1);
        ny=max(2,ceil((varargin{4}-varargin{3})/hmax)+1);
        nz=max(2,ceil((varargin{6}-varargin{5})/hmax)+1);
        end
        x=linspace(varargin{1},varargin{2},nx); 
        y=linspace(varargin{3},varargin{4},ny);
        z=linspace(varargin{5},varargin{6},nz);
        [X,Y,Z]=meshgrid(x,y,z);
    case 4
        if any([length(varargin{1})<=2 ...
                length(varargin{2})<=2 ...
                length(varargin{3})<=2])
            MException('GRID3DPR:WRONGUSE',...
                ['In the 3 argument call the arguments must',...
                ' be vectors defining',...
                ' the coordinates of the points of the bar.']).throw;
        else  % compute h for meshing the base square
            [X,Y,Z]=meshgrid(varargin{1},varargin{2},varargin{3});
            nx=length(varargin{1}); ny=length(varargin{2}); nz=length(varargin{3});
        end
    otherwise
        obj.wrongNumberArguments.throw;
end
P=[reshape(X,1,nx*ny*nz); reshape(Y,1,nx*ny*nz); reshape(Z,1,nx*ny*nz)];
dt=delaunayTriangulation(P'); obj.p=dt.Points'; obj.e=dt.freeBoundary'; obj.t=dt.ConnectivityList';
% tetramesh(dt,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3); pause
obj.t(5,:)=ones(1,size(obj.t,2));  obj=setidsbar(obj); 
end  %%%%%%%%%%% bar %%%%%%%%%%%%%




function ccbar(obj,varargin)  % HU: like bar, but "criss-cross" mesh 
gotnx=0; 
switch nargin
case {7,8,10}
    switch nargin
        case 7; hmax=0.1; case 8; hmax=varargin{7}; 
        case 10; nx=varargin{7}; ny=varargin{8}; nz=varargin{9};gotnx=1; 
    end
    if ~gotnx
    nx=max(2,ceil((varargin{2}-varargin{1})/hmax)+1); ny=max(2,ceil((varargin{4}-varargin{3})/hmax)+1);
    nz=max(2,ceil((varargin{6}-varargin{5})/hmax)+1);
    end
    x=linspace(varargin{1},varargin{2},nx);  y=linspace(varargin{3},varargin{4},ny);
    z=linspace(varargin{5},varargin{6},nz);
    [X,Y,Z]=meshgrid(x,y,z);
otherwise
    obj.wrongNumberArguments.throw;
end
P=[reshape(X,1,nx*ny*nz); reshape(Y,1,nx*ny*nz); reshape(Z,1,nx*ny*nz)];
x1=min(x); x2=max(x); y1=min(y); y2=max(y); z1=min(z); z2=max(z); 
dx=(x2-x1)/(nx-1); dy=(y2-y1)/(ny-1); dz=(z2-z1)/(nz-1); 
for i1=0:nx-2; for i2=0:ny-2; for i3=0:nz-2; % add cube centers and face centers
        xc=x1+dx/2+i1*dx; yc=y1+dy/2+i2*dy; zc=z1+dz/2+i3*dz; 
        cp=[xc;yc;zc]; fxc=[x1+i1*dx; yc; zc]; fyc=[xc; y1+i2*dy; zc]; fzc=[xc; yc; z1+i3*dz]; 
        P=[P cp fxc fyc fzc];
    end; end; end 
for i2=0:ny-2; for i3=0:nz-2; yc=y1+dy/2+i2*dy; zc=z1+dz/2+i3*dz; % add x2 face centers
        fxc=[x2; yc; zc]; P=[P fxc]; end; end;
for i1=0:nx-2; for i3=0:nz-2; xc=x1+dx/2+i1*dx; zc=z1+dz/2+i3*dz; % add y2 face centers
        fyc=[xc; y2; zc]; P=[P fyc]; end; end 
for i1=0:nx-2; for i2=0:ny-2; xc=x1+dx/2+i1*dx; yc=y1+dy/2+i2*dy; % add z2 face centers 
        fzc=[xc; yc; z2]; P=[P fzc]; end; end 
P=P'; P=unique(P,'rows','stable'); P=P'; 
dt=delaunayTriangulation(P'); obj.p=dt.Points'; obj.e=dt.freeBoundary'; obj.t=dt.ConnectivityList';
obj.t(5,:)=ones(1,size(obj.t,2)); obj=setidsbar(obj); 
%tetramesh(dt,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3); pause 
end  %%%%% ccbar %%%%%%%%%%%%%%%%
function unitCube(obj,hmax)
    switch nargin
        case 1 
           hmax=0.2;
        case 2
            % ok
            if ~(isa(hmax,'double')&&isscalar(hmax))
                obj.wrongFormat.throw;
            end
        otherwise
            obj.wrongNumberInputs.throw;
    end
    bar(obj,0,1,0,1,0,1,hmax);      
end        

function ficheraCube(obj)
    % Fichera's cube - test geometry
    % Ficheara's cube is a unit cube with an cutaway of an 1/8
    % volume cube. The number of boundary segements is nine.
    % (c) 2014 by Uwe Prüfert
    obj.p=[    0         0         0
                 0    0.5000         0
                 0    1.0000         0
            0.5000         0         0
            0.5000    0.5000         0
            0.5000    1.0000         0
            1.0000         0         0
            1.0000    0.5000         0
            1.0000    1.0000         0
                 0         0    0.5000
                 0    0.5000    0.5000
                 0    1.0000    0.5000
            0.5000         0    0.5000
            0.5000    0.5000    0.5000
            0.5000    1.0000    0.5000
            1.0000         0    0.5000
            1.0000    0.5000    0.5000
            1.0000    1.0000    0.5000
                 0    0.5000    1.0000
                 0    1.0000    1.0000
            0.5000         0    1.0000
            0.5000    0.5000    1.0000
            0.5000    1.0000    1.0000
            1.0000         0    1.0000
            1.0000    0.5000    1.0000
            1.0000    1.0000    1.0000]';

    obj.e=[    1     2     4     0     1
                 1     4    10     0     3
                 1    10     2     0     2
                 2     3     5     0     1
                 2     5     4     0     1
                 2    10    11     0     2
                 2    11     3     0     2
                 3     6     5     0     1
                 3    11    12     0     2
                 3    12     6     0     5
                 4     5     7     0     1
                 4     7    13     0     3
                 4    13    10     0     3
                 5     6     8     0     1
                 5     8     7     0     1
                 6     9     8     0     1
                 6    12    15     0     5
                 6    15     9     0     5
                 7     8    16     0     4
                 7    16    13     0     3
                 8     9    17     0     4
                 8    17    16     0     4
                 9    15    18     0     5
                 9    18    17     0     4
                10    13    11     0     7
                11    13    14     0     7
                11    14    19     0     8
                11    19    12     0     2
                12    19    20     0     2
                12    20    15     0     5
                13    16    21     0     3
                14    22    19     0     8
                15    20    23     0     5
                15    23    18     0     5
                16    17    24     0     4
                16    24    21     0     3
                17    18    25     0     4
                17    25    24     0     4
                18    23    26     0     5
                18    26    25     0     4
                14    21    22     0     9
                13    21    14     0     9
                19    22    20     0     6
                20    22    23     0     6
                21    24    22     0     6
                22    24    25     0     6
                22    25    23     0     6
                23    25    26     0     6]';

    obj.t=[    6     9    15     8     1
                 5     7    14    13     1
                14     8    16     7     1
                17    24    22    16     1
                24    21    22    16     1
                16    21    14    13     1
                23    20    15    22     1
                 5     4     7    13     1
                23    17    22    15     1
                17     9     8    15     1
                18     9    17    15     1
                23    18    17    15     1
                19    14    20    22     1
                15     8    17    14     1
                15    17    22    14     1
                 3    12     5     6     1
                 5     8    14     7     1
                 3    12    11     5     1
                 2     3    11     5     1
                15    20    14    22     1
                23    25    17    18     1
                23    26    25    18     1
                23    25    22    17     1
                25    24    22    17     1
                16    21    22    14     1
                14     8    17    16     1
                17    16    22    14     1
                14     7    16    13     1
                 4     5    11    13     1
                10     4    11    13     1
                10     2    11     4     1
                 2     5    11     4     1
                10     1     2     4     1
                 5    14    11    13     1
                 6     8    14     5     1
                 6    12     5    14     1
                 5    12    11    14     1
                19    11    12    14     1
                19    12    20    14     1
                15    20    12    14     1
                 6    12    14    15     1
                 6     8    15    14     1]';          
end        

function unitBall(obj,hmax) 
    % unitBall - method that meshes a unit ball
    % obj.unitBall
    % obj.unitBall(hmax)
    % If no hmax given, hmax=0.2 is used 
    % Nate that the used hmax is min(hmax,0.2)
    % to prevent deformed balls.
    % (c) 2014 by Uwe Prüfert.

    switch nargin
        case 1 
            h=0.2;
        case 2
            %
            h=min(0.2,abs(hmax));
        otherwise
            obj.worngNumberInputs.throw
    end
    box=[ -1.0, -1.0, -1.0; ...
           1.0,  1.0,  1.0 ];           

    [obj.p, obj.e, obj.t ]=obj.distmesh(@fd, @fh,h,box,100,[]);

    obj.e=[obj.e;zeros(1,size(obj.e,2));ones(1,size(obj.e,2))];
    obj.t=[obj.t;ones(1,size(obj.t,2))];

    function d=fd( p ) 
        d=sqrt ( sum ( p .^ 2, 2 ) ) - 1.0; 
    end
    function h=fh( p, varargin )                               
        np=size ( p, 1 );
        h=ones ( np, 1 );
    end  
end        

function cylinder(obj,r,l,h,x,y,z)
    % Creates the mesh for a circle geometrie
    % obj.circle() Unitcircle with hmax=0.2
    % obj.circle(R) Circle with Radius=R and hmax=0.2
    % obj.circle(R,hmax)
    % obj.circle(R,hmax,xshift,yshift) 
    % Circle with center=(xshif,yshift)             
    % Code valid for MATLAB R>= 2013a by using new 
    % delaynayTriangulation class. For older Matlab
    % Releases it uses old DelaunayTri Class
    % (c) Uwe Prüfert            

    switch nargin
        case 1
            r=1;
            h=0.2;
            l=1;                    
            x=0;
            y=0;
            z=0;                    
        case 2 
            % r given
            l=1;
            h=r/5;                                   
            x=0;
            y=0;
            z=0;                    
        case 3 
            % r given
            % l given
            h=r/5; 
            x=0;
            y=0;
            z=0;
        case 4
            % r given
            % l given
            % h given
            x=0;
            y=0;
            z=0;
        case 6
            % all given
        otherwise
            obj.wrongNumberInputs.throw;
    end
    L=linspace(0,l,max(2,round(l/h)+1));
    R=linspace(0,r,max(2,round(r/h)+1));
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
    np=size(P,2);
    PZ=[];
    for k=L
        PZ=[PZ,[P;k(ones(1,np))]];
    end


    v=version;
    if str2double(v(1))>= 8
        dt=delaunayTriangulation(PZ');                        
        obj.p=dt.Points';
        obj.e=dt.freeBoundary'; 
        obj.e(5,:)=2*ones(1,size(obj.e,2));
        obj.t=dt.ConnectivityList';  
        obj.t(5,:)=ones(1,size(obj.t,2));  
    else
        dt=DelaunayTri(PZ');                        
        obj.p=dt.X';
        obj.e=dt.freeBoundary'; 
        obj.e(5,:)=2*ones(1,size(obj.e,2));
        obj.t=dt.Triangulation';  
        obj.t(5,:)=ones(1,size(obj.t,2));  
    end
    for k=1:obj.nEdges
        A=obj.p(:,obj.e(1:3,k));
        if A(3,:) == zeros(1,3)
            obj.e(5,k)=1;
        elseif A(3,:) == l(ones(1,3))
            obj.e(5,k)=3;                
        end
    end   
    obj.moveMesh(x,y,z);
end 

function refineMesh(obj,toRefine)
    switch nargin
        case 1
            % uniform refinement
            obj.redRefine;
        case 2
            % not jet implemented, later call the special code
            return
        otherwise
            obj.wrongNumberInputs.throw
    end                        
end
function refineMesh3(obj,toRefine)
    switch nargin
        case 1
            % uniform refinement
            obj.redRefine;
        case 2
            % not jet implemented, later call the special code
            return
        otherwise
            obj.wrongNumberInputs.throw
    end                        
end


function plot(obj,varargin)
    % PLOT method. 
    % obj.plot([patchArguments])
    % Optinal patchArguments may all arguments that can delivered 
    % to patch objects. See also help patch.
    obj.plotFaces([],'FaceColor','none',varargin{:})
    view(3);
    axis equal
end  

function val=evalAtPoint(obj,y,x)
    % evalAtPoint - evaluates y at point x
    %  val=evalAtPoint(obj,y,x)
    val=plotUtils3D().p3interpol(x(1),x(2),x(3),y(:,end),obj.p(1,:),obj.p(2,:),obj.p(3,:));          

    if ~isPointInDomain(obj,x)
        val=NaN;
    end
end   

function [cval,aval,fval]=aCoefficientsMpt(obj,cc,aa,ff)
    %computes the value of the coefficients in the center of every triangle
    % 'symbolic' variables x, y, z  are neccesary for evaluation of string 
    % objects like c='sin(x)' etc.
    np=obj.nPoints; nt=obj.nElements;
    switch class(cc)
        case 'double' % five cases: (i) scalar
            % (ii.1) 3 x 3 diagonal matrix, (ii.2) 3 x 3 full matrix        
            % (iii) vector length np, (iv) vector length nt
            [rows,cols]=size(cc);
            switch max(rows,cols)
                case 1;  cval=cc(ones(3,nt));
                case 3; 
                    if ( norm((ones(3) - diag(3)).*cc,1) == 0)  % (ii.1) 3 x 3 diagonal matrix 
                        c1=cc(1,1); c2=cc(2,2); c3=cc(3,3);
                        cval(1,:)=c1(ones(1,nt)); 
                        cval(2,:)=c2(ones(1,nt)); 
                        cval(3,:)=c3(ones(1,nt));  
                    else  % (ii.2) 3 x 3 full matrix, HU: flipped c! 
                        c1=cc(1,1); c2=cc(1,2); c3=cc(1,3);
                        c4=cc(2,1); c5=cc(2,2); c6=cc(2,3);
                        c7=cc(3,1); c8=cc(3,2); c9=cc(3,3);
    cval(1,:)=c1(ones(1,nt)); cval(2,:)=c2(ones(1,nt)); cval(3,:)=c3(ones(1,nt));   
    cval(4,:)=c4(ones(1,nt)); cval(5,:)=c5(ones(1,nt)); cval(6,:)=c6(ones(1,nt));   
    cval(7,:)=c7(ones(1,nt)); cval(8,:)=c8(ones(1,nt)); cval(9,:)=c9(ones(1,nt));   
                    end
                case nt  % the good one, nothing to do
                    if min(rows,cols)==1; cval=ones(3,1)*cc(:)';
                    elseif (cols==3) || (cols == 9); cval=cc';
                    elseif (rows==3) || (rows == 9); cval=cc;
                    else;  obj.wrongFormat.throw
                    end
                case np % HU, 20-08, c(x,y,z)*diag(1,1,1); 
                    if min(rows,cols)==1; cval=ones(3,1)*obj.point2Center(cc(:));                   
                    else;  obj.wrongFormat.throw
                    end
                otherwise % HU: 3 x (3*np) or 3 x (3*nt) matrix 
      if cols==3*np;  
         cval(1,:)=obj.point2Center(cc(1,1:np)); cval(2,:)=obj.point2Center(cc(1,np+1:2*np));
         cval(3,:)=obj.point2Center(cc(1,2*np+1:3*np));
         cval(4,:)=obj.point2Center(cc(2,1:np)); cval(5,:)=obj.point2Center(cc(2,np+1:2*np));
         cval(6,:)=obj.point2Center(cc(2,2*np+1:3*np));
         cval(7,:)=obj.point2Center(cc(3,1:np)); cval(8,:)=obj.point2Center(cc(3,np+1:2*np));
         cval(9,:)=obj.point2Center(cc(3,2*np+1:3*np));
      else           
         cval(1,:)=cc(1,1:nt); cval(2,:)=cc(1,nt+1:2*nt); cval(3,:)=cc(1,2*nt+1:3*nt);
         cval(4,:)=cc(2,1:nt); cval(5,:)=cc(2,nt+1:2*nt); cval(6,:)=cc(2,2*nt+1:3*nt);
         cval(7,:)=cc(3,1:nt); cval(8,:)=cc(3,nt+1:2*nt); cval(9,:)=cc(3,2*nt+1:3*nt);
      end
                end     
        case 'function_handle'
            % must be a function handle f(x,y,z) of the coefficent function.
            % For evaluating the (string!) coefficient function,
            % we need x,y and z variables "hanging in the air".  
            mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:); y=mpt(2,:); z=mpt(3,:); 
            try cval=feval(cc,x,y,z);
            catch ME; throw(ME); end
        case 'char'
            % must be a single char symbolizing the coefficent function.
            % For evaluating the coefficient function,
            % we need x and y variables "hanging in the air".  
            % Hence, the next warning can be ignored!
            mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:); y=mpt(2,:); z=mpt(3,:); 
            try;  cval=eval(cc);
            catch ME; throw(ME); end
            [rows,cols]=size(cval);
            switch max(rows,cols)
                case 1; cval=cval(ones(3,nt)); % c is a constant like 'pi'                  
                case nt; cval=[cval;cval;cval];
                otherwise;  throw(obj.wrongFormat);
            end
        case 'cell'
            % must be a 3*3 cell array which refers to the entries of a 3*3
            % diffussion matrix. Every entry must be double or a single char symbolizing
            % the part of the coefficent function like in the case 'char'.
            [rows,cols]=size(cc);
            if (rows ~=3) || (cols ~=3); throw(obj.wrongFormat); end
            cval=zeros(9,nt);
            cvalcell=cell(3,3);
            % For evaluating the coefficient function, we need x and y variables "hanging in the air".  
            % Hence, the next warning can be ignored!
            mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:); y=mpt(2,:); z=mpt(3,:); 
            try
                for i=1:3 
                    for j=1 : 3
                        switch class(cc{i,j})
                            case 'char'
                                cvalcell{i,j}=eval(cc{i,j});
                            case 'double'
                                cvalcell{i,j}=cc{i,j};
                            otherwise
                                throw(obj.wrongClass)
                        end
                    end
                end
            catch ME
                throw(ME);
            end
                for i=1:3 
                    for j=1 : 3
                        linIndex=3*(i-1)+j;
                        cvalij=cvalcell{i,j};
                        cvalij=cvalij(:);
                        [rows,cols]=size(cvalij);
                        switch max(rows,cols)
                            case 1 % c is a constant like 'pi'
                                cval(linIndex,:)=cvalij(ones(1,nt));
                            case nt; cval(linIndex,:)=cvalij';
                            otherwise; throw(obj.wrongFormat);
                        end
                    end
                end
        otherwise
            throw(obj.wrongClass)
    end

    % repeat code from case 'c'...
    switch class(aa)
        case 'double'
            [rows,cols]=size(aa);
            switch max(rows,cols)
                case 1
                    aval=aa(ones(1,nt));            
                case np
                    aval=pdeintrp(p,t,aa);
                    aval=[aval;aval];
                case nt
                    % the good one, nothing to do
                    aval=aa;
                otherwise
                    throw(obj.wrongSize)
            end
        case 'function_handle'
           mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:);
            y=mpt(2,:);
            z=mpt(3,:); 
            z= 0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:))); 
            try
                aval=feval(aa,x,y,z);
            catch ME
                throw(ME);
            end
        case 'char'         
            mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:);
            y=mpt(2,:);
            z=mpt(3,:);
            try
                aval=eval(aa);
            catch ME
                throw(ME);
            end
            [rows,cols]=size(aval);
            switch max(rows,cols)
                case 1
                    aval=aval(ones(1,nt));             
                case nt
                    % work already done 
                otherwise
                    throw(finiteElements.wrongFormat);
            end
        case 'cell'
            error('Sorry, Cell array input not jet implemented!')
        otherwise
            throw(obj.wrongClass)
    end 

    switch class(ff)
        case 'double'         
            [rows,cols]=size(ff);
            switch max(rows,cols)
                case 1
                    fval=ff(ones(1,nt));            
                case np
                    fval=pdeintrp(obj.p,obj.t,ff);                 
                case nt
                    % the good one, nothing to do
                    fval=ff;
                otherwise
                    throw(obj.wrongSize)
            end 
        case 'function_handle'
           mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:);
            y=mpt(2,:);
            z=mpt(3,:);

            try
                fval=feval(ff,x,y,z);
            catch ME
                rethrow(ME);
            end

        case 'char'
            mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:);
            y=mpt(2,:);
            z=mpt(3,:);

            try
                fval=eval(ff);
            catch ME
                rethrow(ME);
            end

            [rows,cols]=size(fval);
            switch max(rows,cols)
                case 1
                   fval=fval(ones(1,nt));
                case nt
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

function [bvalvec]=cCoefficients(obj,b)
    %computes the value of the coefficient b in the center of every triangle
    % coefficent can be a vector of dim 2 x 1, or a cell array of dim 2 x 1 

    % 'symbolic' variables x, y  and t are neccesary for evaluation of string 
    % objects like c='sin(x)' etc.



    mpt=1/4*(obj.p(:,obj.t(1,:))...
                    +obj.p(:,obj.t(2,:))...
                    +obj.p(:,obj.t(3,:))...
                    +obj.p(:,obj.t(4,:)));
            x=mpt(1,:);
            y=mpt(2,:);
            z=mpt(3,:);

    % number of points and triangles
    n=obj.nPoints;
    nt=obj.nElements;

    % check the class and  size of b
    if ~((max(size(b))>=3)&&(min(size(b))>=1))
        ME=MException('ccoefficients:wrongCoefficientDefinition',...
            ' b must be a vector');
        throw(ME);
    end
    % two cases:
    % cell-array - entries are strings, or doubles
    if isa(b,'cell')
        b1=b{1};
        b2=b{2};
        b3=b{3};

    elseif isa(b,'double')
        b1=b(1,:);
        b2=b(2,:);
        b3=b(3,:);
    else
        ME=MException('ccoefficients:wrongCoefficientDefinition',...
            ' Wrong coefficient definition');
        throw(ME);
    end


    if isa(b1,'function_handle') || isa(b1,'inline')
        bval=feval(b1,x,y,z);
    elseif isa(b1,'char'),
        bval=eval(b1).*ones(1,n);
    elseif isa(b1,'numeric')
        if length(b1)==n
            % c vektor and defined in p
            bval=b1;
        elseif length(b1)==1,
            % skalar
            bval=b1*ones(nt,1);
        elseif length(b1)==nt,
            bval=b1;
        else
            ME=MException('ccoefficients:wrongSize','wrong sized b(1)');
            throw(ME);
        end
    elseif isa(b1,'inline')
         bval=b(x,y,z);
    else
        ME=MException('ccoefficients:wrongSize','wrong formated b(1)');
        throw(ME);
    end


    if length(bval) == nt
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
        bval=feval(b2,x,y,z);
    elseif isa(b2,'char'),
        bval=eval(b2).*ones(1,n);
    elseif isa(b2,'numeric')
        if length(b2)==n
            % c vektor and defined in p
            bval=b2;
        elseif length(b2)==1,
            % skalar
            bval=b2*ones(nt,1);
        elseif length(b2)==nt,
            bval=b2;
        else
            ME=MException('ccoefficients:wrongSize','wrong sized b(1)');
            throw(ME);
        end
    elseif isa(c,'inline')
         bval=b(x,y,z);
    else
        ME=MException('ccoefficients:wrongSize','wrong formated b(1)');
        throw(ME);
    end

    if length(bval) == nt
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


    if isa(b3,'function_handle') || isa(b3,'inline')
        bval=feval(b3,x,y,z);
    elseif isa(b2,'char'),
        bval=eval(b2).*ones(1,n);
    elseif isa(b2,'numeric')
        if length(b2)==n
            % c vektor and defined in p
            bval=b2;
        elseif length(b2)==1,
            % skalar
            bval=b2*ones(nt,1);
        elseif length(b2)==nt,
            bval=b2;
        else
            ME=MException('ccoefficients:wrongSize','wrong sized b(1)');
            throw(ME);
        end
    elseif isa(b3,'inline')
         bval=b3(x,y,z);
    else
        ME=MException('ccoefficients:wrongSize','wrong formated b(1)');
        throw(ME);
    end




    if length(bval) == nt
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
    bvalvec= bvalvec';
end

function b=isPointInElement(obj,pt)
    % isPointInTetraeder -  gives back a bool vector of length nElements.
    % b ist true if the point pt is i the tetraeder of a given
    % mesh object.
    % b=isPointInTetreder(obj,pt)  
    % Algortihm: Transform the point wrt the transoformation of 
    % the tetraeder into the unit tetraeder  and
    % decide on basis of the relation pt_x < 0 pt_y < 0 ...
    p=obj.p;
    t=obj.t;

    p1=p(:,(t(1,:)));
    p2=p(:,(t(2,:)));
    p3=p(:,(t(3,:)));
    p4=p(:,(t(4,:)));

    x21=p2(1,:)-p1(1,:);
    x31=p3(1,:)-p1(1,:);
    x41=p4(1,:)-p1(1,:);

    y21=p2(2,:)-p1(2,:);
    y31=p3(2,:)-p1(2,:); 
    y41=p4(2,:)-p1(2,:); 

    z21=p2(3,:)-p1(3,:);
    z31=p3(3,:)-p1(3,:); 
    z41=p4(3,:)-p1(3,:); 


    J= ( (x21.*y31-x31.*y21).*z41...
        -(x21.*y41-x41.*y21).*z31...
        +(x31.*y41-x41.*y31).*z21); 

    % transform the point into the unit-triangle and
    % decide
    ptx=pt(1)-p1(1,:);
    pty=pt(2)-p1(2,:);
    ptz=pt(3)-p1(3,:);

    % to compute     the inverse of linear transformation
    xi_x=(y31.*z41-y41.*z31)./J; 
    eta_x=(y41.*z21-y21.*z41)./J;
    zeta_x= (y21.*z31-y31.*z21)./J; % zero for prisms

    xi_y=(x41.*z31-x31.*z41)./J ;
    eta_y=(x21.*z41-x41.*z21)./J;  
    zeta_y=(x31.*z21-x21.*z31)./J; % zero for prisms

    xi_z=(x31.*y41-x41.*y31)./J;
    eta_z=(x41.*y21-x21.*y41)./J;
    zeta_z=(x21.*y31-x31.*y21)./J;



    ptux=(xi_x.*ptx+xi_y.*pty+xi_z.*ptz);
    ptuy=(eta_x.*ptx+eta_y.*pty+eta_z.*ptz); 
    ptuz=(zeta_x.*ptx+zeta_y.*pty+zeta_z.*ptz); 

    b=~(sum([ptux;ptuy;ptuz])>1 | ptux<0 | ptuy<0 | ptuz<0);
end        
end

methods(Static) 
function[qval,gval,hval,rval]=bCoefficients(p,b)
    % compute the boundary coefficients
    % 
    % q,g,h,r are srings defining everything that can be 
    % evaluated by eval. The
    % independent variables must be named by x, y, z
    % (c) 2013 by Uwe Prüfert

    x=sum(p(1,:))/size(p,2); %#ok<*NASGU>
    y=sum(p(2,:))/size(p,2);
    z=sum(p(3,:))/size(p,2);    


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

methods(Access=private)       
function obj=redRefine(obj)
   % redRefine - refines 3D tetraeder mesh uniformly by
   % "red" bisection, i.e. each tetraeder will be refined by
   % adding new knodes on the middle of each edge and add eight 
   % new tetraeders build from them to the mesh.
   % (c) 2014 by Uwe Prüfert. The code bases on an idea taken
   % from the Bachelorthesis of J.F. Kemetmüller.


    % zero-vector to fill the edge matrix (it contains
    % placeholders for the fourth edge point used by prism
    % elements
    n=zeros(1,obj.nEdges);

    % do not call the method too often...
    npts=obj.nPoints;

    % Use the idea of Kemetmüller to compute the conectivity
    % matrices
    [element2edges,edge2nodes,boundary2edges] ...
       =obj.provideEdgeData();

    obj.p=[obj.p 0.5*(obj.p(:,edge2nodes(:,1))+obj.p(:,edge2nodes(:,2)))];


    % The original code was wrong (or not suitable for positive 
    % oriented elements). I change the order in the code for the 
    % 6th and 7th new element to correct (adapt ?) it.   
    obj.t=[[obj.t(1,:);npts+element2edges([1,2,3],:);obj.t(5,:)] ,...
            [npts+element2edges(1,:);obj.t(2,:);npts+element2edges([4,5],:);obj.t(5,:)],...
            [npts+element2edges([2,4],:);obj.t(3,:);npts+element2edges(6,:);obj.t(5,:)],...
            [npts+element2edges([3,5,6],:);obj.t(4,:);obj.t(5,:)],...
            [npts+element2edges([1,2,3,5],:);obj.t(5,:)],...
            [npts+element2edges([1,4,2,5],:);obj.t(5,:)],...
            [npts+element2edges([2,5,4,6],:);obj.t(5,:)],...
            [npts+element2edges([2,3,5,6],:);obj.t(5,:)]]; 

    obj.e=[[obj.e(1,:);npts+boundary2edges([3,2],:);n;obj.e(5,:)],...
            [obj.e(2,:);npts+boundary2edges([1,3],:);n;obj.e(5,:)],...
            [obj.e(3,:);npts+boundary2edges([2,1],:);n;obj.e(5,:)],...
            [npts+boundary2edges([1 2 3],:);n;obj.e(5,:)]];
end
 function obj=redRefine3(obj, sel) % HU 

    n=zeros(1,obj.nEdges);
    npts=obj.nPoints;
    [element2edges,edge2nodes,boundary2edges] ...
       =obj.provideEdgeData();

    obj.p=[obj.p 0.5*(obj.p(sel,edge2nodes(sel,1))+obj.p(sel,edge2nodes(sel,2)))];


    % The original code was wrong (or not suitable for positive 
    % oriented elements). I change the order in the code for the 
    % 6th and 7th new element to correct (adapt ?) it.   
    obj.t=[[obj.t(1,:);npts+element2edges([1,2,3],:);obj.t(5,:)] ,...
            [npts+element2edges(1,:);obj.t(2,:);npts+element2edges([4,5],:);obj.t(5,:)],...
            [npts+element2edges([2,4],:);obj.t(3,:);npts+element2edges(6,:);obj.t(5,:)],...
            [npts+element2edges([3,5,6],:);obj.t(4,:);obj.t(5,:)],...
            [npts+element2edges([1,2,3,5],:);obj.t(5,:)],...
            [npts+element2edges([1,4,2,5],:);obj.t(5,:)],...
            [npts+element2edges([2,5,4,6],:);obj.t(5,:)],...
            [npts+element2edges([2,3,5,6],:);obj.t(5,:)]]; 

    obj.e=[[obj.e(1,:);npts+boundary2edges([3,2],:);n;obj.e(5,:)],...
            [obj.e(2,:);npts+boundary2edges([1,3],:);n;obj.e(5,:)],...
            [obj.e(3,:);npts+boundary2edges([2,1],:);n;obj.e(5,:)],...
            [npts+boundary2edges([1 2 3],:);n;obj.e(5,:)]];
end

function [element2edges,edge2nodes,boundary2edges]=provideEdgeData(obj)
    %  provideEdgeData - Computes the conectivity matrices for
    %  refine tetraeders and identify the "right" points on
    %  the new edges end tetraeders. 
    % (c) 2014 by Uwe Prüfert. The code is an adaption of
    % Kemetmüller's code.

    elements=obj.t(1:4,:)';
    boundary=obj.e(1:3,:)';

    bEInd=cumsum([obj.nElements*6,3*obj.nEdges]);

    edgeorder=[1,2;...
                 1,3;...
                 1,4;...
                 2,3;...
                 2,4;...
                 3,4];

    boundaryedgeorder=[2,3;...
                         1,3;...
                         1,2];                             

    edges=sort(reshape(elements(:,edgeorder),6*obj.nElements,2),2);

    edges(bEInd(1)+1:bEInd(2),:)=...
        sort(reshape(boundary(:,boundaryedgeorder),[],2),2);

    [edge2nodes,~,J]=unique(edges,'rows');
    element2edges=reshape(J(1:obj.nElements*6),obj.nElements,6)';            

    boundary2edges=reshape(J(bEInd(1)+1:bEInd(2)),obj.nEdges,3)';
end
end

methods(Static,Access=private,Hidden) 
function [p,e,t]=distmesh(fd, fh, h0, box, iteration_max, pfix, varargin )
    % DISTMESH_3D 3D Mesh Generator using Distance Functions.
    % Adapted for OOPDE by Uwe Prüfert.
    %  Copyright:
    %
    %    (C) 2004 Per-Olof Persson. 
    %    See COPYRIGHT.TXT for details.
    %
    %   
    %  Parameters:
    %
    %    Input, function_handle FD  signed distance function d(x,y).
    %           function_handle FH  scaled edge length function h(x,y).
    %           double              H0, the initial smallest edge length.
    %           double              BOX(3,2), a bounding box for the region.
    %           double              ITERATION_MAX, the maximum number of iterations.
    %                               The iteration might terminate sooner than this
    %                               limit, if the program decides that the mesh 
    %                               has converged.
    %           double PFIX(NFIX,3) the coordinates of nodes that are
    %                               required to be included in the mesh.
    %           VARARGIN,           additional parameters that can be passed to FD
    %
    %    Output, double P(N,3)      node coordinates.
    %            double T(NT,4)     tetrahedron indices.
    ptol=0.001; 
    ttol=0.1; 
    L0mult=1.1; 
    deltat=0.1; 

    geps=0.1*h0;
    deps=sqrt(eps)*h0;
    iteration=0;
    iteration_max=max(iteration_max,1);

    cbox=cell(1,3);
    for k=1 : 3
        cbox{k}=box(1,k):h0:box(2,k);
    end
    pp=cell(1,3);
    [pp{:}]=ndgrid(cbox{:});
%             p=zeros(prod(size(pp{1})),3);
    p=zeros(numel(pp{1}),3);
    for k=1 : 3
        p(:,k)=pp{k}(:);
    end
    %
    % 2. Remove points outside the region, apply the rejection method.
    %
    p=p (feval(fd,p,varargin{:}) < geps,:);
    r0=feval (fh,p);
    p=p(rand(size(p,1),1)<min(r0)^3./r0.^3, : );
    p=[pfix;p];  

    N=size ( p, 1 );

    %   count=0;
    p0=inf;

    while ( iteration < iteration_max )  
        iteration=iteration + 1;
        if ( ttol * h0 < max (sqrt(sum((p-p0).^2,2)) ) )
            p0=p;
            t=delaunayn ( p );

            pmid=zeros(size(t,1),3);
            for k=1 : 3+1
                pmid=pmid + p(t(:,k),:) / (3+1);
            end
            %
            %  Only keep those simplices in the new triangulation for which
            %  the centroids is inside the region, or not too far outside.
            %
            t=t ( feval ( fd, pmid, varargin{:} ) < -geps, : );
            %
            %  4. Describe each edge by a unique pair of nodes.
            %
            pair=zeros(0,2);
            localpairs=nchoosek(1:3+1,2);
            for k=1 : size(localpairs,1)
                pair=[pair;t(:,localpairs(k,:))];
            end
            pair=unique(sort(pair,2),'rows'); 
        end

        bars=p(pair(:,1),:)-p(pair(:,2),:);
        L=sqrt(sum(bars.^2,2));
        L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
        L0=L0*L0mult*(sum(L.^3)/sum(L0.^3))^(1/3);
        F=max(L0-L,0);
        Fbar=[bars,-bars].*repmat(F./L,1,2*3);
        dp=full(sparse(pair(:,[ones(1,3),2*ones(1,3)]), ...
                ones(size(pair,1),1)*[1:3,1:3], ...
                Fbar,N,3));
        dp(1:size(pfix,1),:)=0;
        p=p + deltat * dp;

        for steps=1 : 2

            d=feval ( fd, p, varargin{:} );
            indx=( 0 < d );
            gradd=zeros ( sum(indx), 3 );

            for k=1 : 3
                a=zeros(1,3);
                a(k)=deps;
                d1x=feval ( fd, p(indx,:)+ones(sum(indx),1)*a, varargin{:} );
                gradd(:,k)=( d1x - d(indx) ) / deps;
            end

            p(indx,:)=p(indx,:) - d(indx) * ones(1,3) .* gradd;

        end

        maxdp=max ( deltat * sqrt ( sum(dp(d < -geps, : ).^2, 2 ) ) );

        if ( maxdp < ptol * h0 )
            break;
        end

    end
    e=surf2tri ( p, t );
    p=p';
    e=e';
    t=t';

    function tri=surf2tri ( p, t )                
        %
        faces=[  t(:,[1,2,3]);
                   t(:,[1,2,4]);
                   t(:,[1,3,4]);
                   t(:,[2,3,4])];
        node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
        faces=sort(faces,2);
        [~,iindx,jindx]=unique(faces,'rows');
        vec=histc(jindx,1:max(jindx));
        qx=find(vec==1);
        tri=faces(iindx(qx),:);
        node4=node4(iindx(qx));
        %
        % Orientation
        %
        v1=p(tri(:,2),:)-p(tri(:,1),:);
        v2=p(tri(:,3),:)-p(tri(:,1),:);
        v3=p(node4,:)-p(tri(:,1),:);
        iindx=find(dot(cross(v1,v2,2),v3,2)>0);
        tri(iindx,[2,3])=tri(iindx,[3,2]);
    end
end

function d=ddiff (d1,d2)
    % DDIFF returns the signed distance to a region that is the difference of two regions.
    % Taken from
    %    Per-Olof Persson and Gilbert Strang,
    %    A Simple Mesh Generator in MATLAB,
    %    SIAM Review,
    %    Volume 46, Number 2, June 2004, pages 329-345.
    %
    %  Parameters:
    %
    %    Input, real D1, D2, the signed distances to region 1 and 2.
    %
    %    Output, real D, the signed distance to the region formed by
    %    removing from region 1 its intersection with region 2.
    %
    d=max ( d1, -d2 );             
end

function d=dintersect ( d1, d2 )
    % DINTERSECT sets the signed distance to the intersection of two regions.
    % Taken from
    %    Per-Olof Persson and Gilbert Strang,
    %    A Simple Mesh Generator in MATLAB,
    %    SIAM Review,
    %    Volume 46, Number 2, June 2004, pages 329-345.
    %
    %  Parameters:
    %
    %    Input, real D1, D2, the signed distance to each of the regions.
    %
    %    Output, real D, the signed distance to the region formed by the
    %    intersection of the two regions.
    %
    d=max(d1,d2);         
end
end  
end

