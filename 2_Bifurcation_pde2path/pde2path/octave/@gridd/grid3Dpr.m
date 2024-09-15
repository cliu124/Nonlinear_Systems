classdef grid3Dpr < gridd & plotUtils3D
    %grid3D 3D grid for used with FEM-OO Toolbox
    % Creates prism element meshes. 
    % (c) 2013 Uwe Prüfert
    properties(Access = public)
        basis;
    end
    methods(Access = public)        
        function obj = grid3Dpr()
            % constructor method
            % only empty call allowed
            % obj = grid3Dpr()
            
            % only  in geometry
            obj.nPointsInElements = 6;
            switch nargin
                case 0
                    % empty object
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end 
        
        % some geometries
        function unitCube(obj,h)
            switch nargin
                case 1
                    h = 0.1;
                case 2   
                    %
                otherwise
            end
            base = grid2D();
            base.unitSquare(h);
            obj.extrude(base, linspace(0,1,ceil(2/h)+1));
            
            % Correct the Boundary numbering...
            % because the unitcube has its own system
            for k = 1:obj.nEdges
                A = obj.p(:,obj.e(1:3,k));
                if A(1,:) == zeros(1,3)
                    obj.e(5,k) = 2;
                elseif A(1,:) == ones(1,3)
                    obj.e(5,k) = 4;
                elseif A(2,:) == zeros(1,3)
                    obj.e(5,k) = 3;
                elseif A(2,:) == ones(1,3)
                    obj.e(5,k) = 5;
                elseif A(3,:) == zeros(1,3)
                    obj.e(5,k) = 1;
                elseif A(3,:) == ones(1,3)
                    obj.e(5,k) = 6;                
                end
            end      

        end
        
        function bar(obj,varargin)
            % method that meshes a 3D bar.
            % obj.bar(xmin,xmax,ymin,ymax,zmin,zmax)
            % obj.bar(xmin,xmax,ymin,ymax,zmin,zmax,hmax)
            % obj.bar(xvec,yvec,zvec)
            % xvec, ybvec zvec must be vectors with length >2.
            % This call is intentionally for using different meshwides in
            % x-y and z direction.
            % (c) 2013 by Uwe Prüfert 
            
            base = grid2D();
            switch nargin
                case {4}                   
                    if any([...
                            length(varargin{1})<=2 ...
                            length(varargin{2})<=2 ...
                            length(varargin{3})<=2])
                        MException('GRID3DPR:WRONGUSE',...
                            ['In the 3 argument call the arguments must',...
                            ' be vectors defining',...
                            ' the coordinates of the points of the bar.']).throw;
                    else
                        % compute h for meshing the base square
                        h = min([varargin{1}(2:end)-varargin{1}(1:end-1),...
                            varargin{2}(2:end)-varargin{2}(1:end-1)]);
                    end
                    
                    base.square(varargin{1}(1),varargin{1}(end),...
                        varargin{2}(1),varargin{2}(end),h);                    
                    obj.extrude(base,varargin{3});
                case {7,8}
                    switch nargin
                        case 8
                            h = varargin{7};
                        case 7
                            h = 0.1;                 
                        otherwise
                            obj.wrongNumberInputs.throw
                    end
                   
                    base.square(varargin{1},varargin{2},varargin{3},varargin{4},h);
                    obj.extrude(base,linspace(varargin{5},varargin{6},...
                        ceil((varargin{6}-varargin{5})/h)+1));
                otherwise
                    MException('GRID3DPR:WRONGNUMBEROFARGUMENTS',...
                        [' The number of arguments must be 3, 6 or 7.',...
                          'See help grid3Dpr.bar.']).throw;
            end
            
            

        end
        
        function cylinder(obj,r,l,rmax,lmax)
            switch nargin
                case 1
                    r = 1;
                    l = 1;
                    rmax = 0.2;
                    lmax = 0.3;                 
                case 3
                    rmax = 0.2;
                    lmax = 0.3;
                case 4
                    lmax = 0.7*rmax;
                case 5
                    % full input list...
                otherwise
                    obj.wrongNumberInputs.throw;
            end
            obj.basis = grid2D();
            obj.basis.circle(r,rmax); 
            obj.extrude(obj.basis,linspace(0,l,ceil(l/lmax)+1));
        end
        
        function rail(obj,l,lmax)
            switch nargin
                case 1
                    l = 3;
                     lmax = 0.3;   
                case 2
                    lmax = l/10;
                case 3                     
                    % full input list...
                otherwise
                    obj.wrongNumberInputs.throw;
            end
            base = grid2D();
            base.doubleT();
            base.refinemesh();
            base.refinemesh(); 
            problem.grid.extrude(base,linspace(0,l,ceil(l/lmax)+1));
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
              
        function [cval,aval,fval] = aCoefficientsMpt(obj,cc,aa,ff)
            %computes the value of the coefficients in the center of every triangle
            % 'symbolic' variables x, y  are neccesary for evaluation of string 
            % objects like c = 'sin(x)' etc.
            % Restrictions
            % NOT jet implemented:
            % changeLog:
            % 03/03/2015
            % cell array input added
            % f ull matrix c input added


            % p = obj.p;
            % t = obj.t;
            np = obj.nPoints;
            nt = obj.nElements;

            switch class(cc)
%                 case 'function_handle'
%                     mpt = 1/3*(obj.p(:,obj.t(1,:))...
%                             +obj.p(:,obj.t(2,:))...
%                             +obj.p(:,obj.t(3,:)));
%                     x = mpt(1,:);
%                     y = mpt(2,:);
%                     z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
%                     try
%                         c = feval(cc,x,y,z);
%                     catch ME
%                         rethrow(ME);
%                     end  
%                     cval(1,:) = c(1,1:nt)); 
%                     cval(2,:) = c2(ones(1,nt)); 
%                     cval(3,:) = c3(ones(1,nt));  
                    
                case 'double'
                    % five case
                    % (i) scalar
                    % (ii.1) 3 x 3 diagonal matrix     
                    % (ii.2) 3 x 3 full matrix        
                    % (iii) vector length np
                    % (iv) vector length nt
                    [rows,cols] = size(cc);
                    switch max(rows,cols)
                        case 1
                            cval = cc(ones(3,nt));
                        case 3
                           if ( norm((ones(3) - diag(3)).*cc,1) == 0)
                               % (ii.1) 3 x 3 diagonal matrix 
                                c1 = cc(1,1); c2 = cc(2,2); c3 = cc(3,3);
                                cval(1,:) = c1(ones(1,nt)); 
                                cval(2,:) = c2(ones(1,nt)); 
                                cval(3,:) = c3(ones(1,nt));   
                           else 
                                % (ii.2) 3 x 3 full matrix   
                                c1 = cc(1,1); c2 = cc(2,1); c3 = cc(3,1);
                                c4 = cc(1,2); c5 = cc(2,2); c6 = cc(3,2);
                                c7 = cc(1,3); c8 = cc(2,3); c9 = cc(3,3);
                                cval(1,:) = c1(ones(1,nt)); 
                                cval(2,:) = c2(ones(1,nt)); 
                                cval(3,:) = c3(ones(1,nt));   
                                cval(4,:) = c4(ones(1,nt)); 
                                cval(5,:) = c5(ones(1,nt)); 
                                cval(6,:) = c6(ones(1,nt));   
                                cval(7,:) = c7(ones(1,nt)); 
                                cval(8,:) = c8(ones(1,nt)); 
                                cval(9,:) = c9(ones(1,nt));   
                           end
                        case nt
                            % the good one, nothing to do
                            if min(rows,cols)==1
                                cval = ones(3,1)*cc(:)';
                            elseif (cols==3) || (cols == 9)
                                cval = cc';
                            elseif (rows==3) || (rows == 9)
                                cval = cc;
                            else
                                obj.wrongFormat.throw
                            end
                                
                        case np
                            
                            MException('grid3Dpr:aCoefficients',...
                                ['Not jet implemented, try to use a function',... 
                                ' for the midpoints of elements']).throw
                            
                        otherwise
                             obj.wrongFormat.throw
                    end
                case 'char'
                    % must be a single char symbolizing
                    % the coefficent function.
                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".  
                    % Hence, the next warning can be ignored!
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        cval = eval(cc);
                    catch ME
                        throw(ME);
                    end
                    [rows,cols] = size(cval);
                    switch max(rows,cols)
                        case 1
                            % c is a constant like 'pi'
                            cval = cval(ones(3,nt));
                        case nt
                            cval = [cval;cval;cval];
                        otherwise
                            throw(obj.wrongFormat);
                    end
                case 'cell'
                    % must be a 3*3 cell array which refers to the entries of a 3*3
                    % diffussion matrix. Every Entry must be double or a single char symbolizing
                    % the part of the coefficent function like in the case
                    % 'char'.
                    
                    [rows,cols] = size(cc);
                    if (rows ~=3) || (cols ~=3)
                        throw(obj.wrongFormat);
                    end
                    cval = zeros(9,nt);
                    cvalcell = cell(3,3);
                    
                    
                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".  
                    % Hence, the next warning can be ignored!
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        for i = 1:3 
                            for j = 1 : 3
                                switch class(cc{i,j})
                                    case 'char'
                                        cvalcell{i,j} = eval(cc{i,j});
                                    case 'double'
                                        cvalcell{i,j} = cc{i,j};
                                    otherwise
                                        throw(obj.wrongClass)
                                end
                            end
                        end
                    catch ME
                        throw(ME);
                    end
                        for i = 1:3 
                            for j = 1 : 3
                                linIndex = 3*(i-1)+j;
                                cvalij = cvalcell{i,j};
                                cvalij = cvalij(:);
                                [rows,cols] = size(cvalij);
                                switch max(rows,cols)
                                    case 1
                                        % c is a constant like 'pi'
                                        cval(linIndex,:) = cvalij(ones(1,nt));
                                    case nt
                                        cval(linIndex,:) = cvalij';
                                    otherwise
                                        throw(obj.wrongFormat);
                                end
                            end
                        end
                otherwise
                    throw(obj.wrongClass)
            end
            % repeat code from case 'c'...
            switch class(aa)
                case 'double'
                    [rows,cols] = size(aa);
                    switch max(rows,cols)
                        case 1
                            aval = aa(ones(1,nt));            
                        case np
                            aval = pdeintrp(p,t,aa);
                            aval = [aval;aval];
                        case nt
                            % the good one, nothing to do
                            aval = aa;
                        otherwise
                            throw(obj.wrongSize)
                    end
                case 'char'         
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:))); 
                    try
                        aval = eval(aa);
                    catch ME
                        throw(ME);
                    end
                    [rows,cols] = size(aval);
                    switch max(rows,cols)
                        case 1
                            aval = aval(ones(1,nt));             
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
                    [rows,cols] = size(ff);
                    switch max(rows,cols)
                        case 1
                            fval = ff(ones(1,nt));            
                        case np
                            fval = pdeintrp(obj.p,obj.t,ff);                 
                        case nt
                            % the good one, nothing to do
                            fval = ff;
                        otherwise
                            throw(obj.wrongSize)
                    end
                case 'function_handle'
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        fval = feval(ff,x,y,z);
                    catch ME
                        rethrow(ME);
                    end                     
                case 'char'
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));

                   try
                        fval = eval(ff);
                    catch ME
                        rethrow(ME);
                    end

                    [rows,cols] = size(fval);
                    switch max(rows,cols)
                        case 1
                           fval = fval(ones(1,nt));
                        case nt
                            % work already done
                        otherwise
                            throw(obj.wrongFormat);
                    end
                case 'cell'
                    error('Sorry, Cell array input makes here no sense.')
                otherwise
                    throw(obj.wrongClass)
            end
        end
  
        function [bvalvec] = cCoefficients(g,b)
            %computes the value of the coefficient b in the center of every triangle
            % coefficent can be a vector of dim 2 x 1, or a cell array of dim 2 x 1 

            % 'symbolic' variables x, y  and t are neccesary for evaluation of string 
            % objects like c = 'sin(x)' etc.
            p = g.p;
            t = g.t;

            x = p(1,:);
            y = p(2,:);
            z = p(3,:);

            % number of points and triangles
            n = g.nPoints;
            nt = g.nElements;

            % check the class and  size of b
            if ~((max(size(b))>=2)&&(min(size(b))>=1))
                ME = MException('ccoefficients:wrongCoefficientDefinition',...
                    ' b must be a vector');
                throw(ME);
            end
            % two cases:
            % cell-array - entries are strings, or doubles
            if isa(b,'cell')
                b1 = b{1};
                b2 = b{2};              
                b3 = b{3};
            elseif isa(b,'double')
                b1 = b(1,:);
                b2 = b(2,:); 
                b3 = b(3,:);
            else
                ME = MException('ccoefficients:wrongCoefficientDefinition',...
                    ' Wrong coefficient definition');
                throw(ME);
            end


            if isa(b1,'function_handle') || isa(b1,'inline')
                bval = feval(b1,x,y);
            elseif isa(b1,'char'),
                bval = eval(b1).*ones(1,n);
            elseif isa(b1,'numeric')
                if length(b1)==n
                    % c vektor and defined in p
                    bval = b1;
                elseif length(b1)==1,
                    % skalar
                    bval = b1*ones(nt,1);
                elseif length(b1)==nt,
                    bval = b1;
                else
                    ME = MException('ccoefficients:wrongSize','wrong sized b(1)');
                    throw(ME);
                end
            elseif isa(b1,'inline')
                 bval = b(x,y);
            else
                ME = MException('ccoefficients:wrongSize','wrong formated b(1)');
                throw(ME);
            end
            if length(bval) == nt
                % b(1) is a vektor and defined in center of mass of triagle
            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = g.point2Center(bval);  
            end
            dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
            % first column of bvalvec
            bvalvec = bval;

            if isa(b2,'function_handle') || isa(b2,'inline')
                bval = feval(b2,x,y);
            elseif isa(b2,'char'),
                bval = eval(b2).*ones(1,n);
            elseif isa(b2,'numeric')
                if length(b2)==n
                    % c vektor and defined in p
                    bval = b2;
                elseif length(b2)==1,
                    % skalar
                    bval = b2*ones(nt,1);
                elseif length(b2)==nt,
                    bval = b2;
                else
                    ME = MException('ccoefficients:wrongSize','wrong sized b(1)');
                    throw(ME);
                end
            elseif isa(b2,'inline')
                 bval = b(x,y);
            else
                ME = MException('ccoefficients:wrongSize','wrong formated b(1)');
                throw(ME);
            end

            if length(bval) == nt
                % b(1) is a vektor and defined in center of mass of triagle

            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = g.point2Center(bval);  
            end
            dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
            bvalvec = [bvalvec,bval];
            
            if isa(b3,'function_handle') || isa(b2,'inline')
                bval = feval(b2,x,y);
            elseif isa(b3,'char'),
                bval = eval(b2).*ones(1,n);
            elseif isa(b3,'numeric')
                if length(b3)==n
                    % c vektor and defined in p
                    bval = b3;
                elseif length(b3)==1,
                    % skalar
                    bval = b3*ones(nt,1);
                elseif length(b3)==nt,
                    bval = b3;
                else
                    ME = MException('ccoefficients:wrongSize','wrong sized b(1)');
                    throw(ME);
                end
            elseif isa(b3,'inline')
                 bval = b(x,y);
            else
                ME = MException('ccoefficients:wrongSize','wrong formated b(1)');
                throw(ME);
            end
            if length(bval) == nt
                % b(1) is a vektor and defined in center of mass of triagle

            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = g.point2Center(bval);  
            end

            dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
            bvalvec = [bvalvec,bval]';
            
        end 
          
        function b = isPointInElement(obj,pt)
            % isPointInTetraeder -  gives back a bool vector of length nElements.
            % b ist true if the point pt is i the tetraeder of a given
            % mesh object.
            % b = isPointInTetreder(obj,pt)  
            % Algortihm: Transform the point wrt the transoformation of 
            % the tetraeder into the unit tetraeder  and
            % decide on basis of the relation pt_x < 0 pt_y < 0 ...
            p = obj.p;
            t = obj.t;
            
            p1 = p(:,(t(1,:)));
            p2 = p(:,(t(2,:)));
            p3 = p(:,(t(3,:)));
            p4 = p(:,(t(4,:)));
            
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            x41 = p4(1,:)-p1(1,:);
            
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:); 
            y41 = p4(2,:)-p1(2,:); 
            
            z21 = p2(3,:)-p1(3,:);
            z31 = p3(3,:)-p1(3,:); 
            z41 = p4(3,:)-p1(3,:); 
            
            
            J =  ( (x21.*y31-x31.*y21).*z41...
                -(x21.*y41-x41.*y21).*z31...
                +(x31.*y41-x41.*y31).*z21); 
            
            % transform the point into the unit-triangle and
            % decide
            ptx = pt(1)-p1(1,:);
            pty = pt(2)-p1(2,:);
            ptz = pt(3)-p1(3,:);
            
            % to compute     the inverse of linear transformation
            xi_x = (y31.*z41-y41.*z31)./J; 
            eta_x = (y41.*z21-y21.*z41)./J;
            zeta_x =  (y21.*z31-y31.*z21)./J; % zero for prisms
            
            xi_y = (x41.*z31-x31.*z41)./J ;
            eta_y = (x21.*z41-x41.*z21)./J;  
            zeta_y = (x31.*z21-x21.*z31)./J; % zero for prisms
            
            xi_z = (x31.*y41-x41.*y31)./J;
            eta_z = (x41.*y21-x21.*y41)./J;
            zeta_z = (x21.*y31-x31.*y21)./J;
 
            ptx = pt(1)-p1(1,:);
            pty = pt(2)-p1(2,:);
            ptz = pt(3)-p1(3,:);
            ptux = (xi_x.*ptx+xi_y.*pty+xi_z.*ptz);
            ptuy = (eta_x.*ptx+eta_y.*pty+eta_z.*ptz);
            ptuz = (zeta_z.*ptz);
            
            b = ~(sum([ptux;ptuy])>1 | ptuz>1 | ptux<0 | ptuy<0 | ptuz<0);
        end        
    end
    
    
    methods(Access = protected)
        function extrude(obj,grid2d,z)
        %extrude(obj,grid2d,z)
            %  Extrudes a 2D geometry in z direction. grid2d must be a valid GRIDD2D
            %  object defining a 2D geometry object. The result is a
            %  6-NODE WEDGE ELEMENT (aka prism element) mesh.
            if ~isa(grid2d,'grid2D')|| min(size(grid2d.p)) ~= 2
                throw(obj.wrongClass);
            end
            nSlices = length(z);
            nPts = grid2d.nPoints;
            nTri =  grid2d.nElements;
            if num2str(nTri*(nSlices-1))>1e6
                fprintf(['Warning: A mesh with ',num2str(nTri*(nSlices-1)),...
                    ' elements and ',num2str(nSlices*nPts),...
                    ' points will now be generated.\n']);
            end
            obj.p = zeros(3,nPts*nSlices);
            block = expandFirstSlice(grid2d.t,nPts);
            obj.t = block;
            
               % the edges
            % initialize with triangel data from 2D object
            obj.e = grid2d.t(1:3,:);
            % base gets the subdomain numbers from 2D object
            obj.e(5,1:nTri) = grid2d.t(4,:) ;
            nBdConds = max(obj.e(5,:));
            EdgeBlock = expandBoundary(grid2d.e,nPts,nBdConds);
            obj.e = [obj.e EdgeBlock];            
             
            % wedge and points-loop
            for k = 1:nSlices-2
                block = block+nPts;
                EdgeBlock(1:4,:) = EdgeBlock(1:4,:)+nPts;
                obj.t = [obj.t block];
                obj.p(:,nPts*(k-1)+1:nPts*k) = ...
                    [grid2d.p;z(k)*ones(1,nPts)];
                obj.e = [obj.e [EdgeBlock(1:4,:);EdgeBlock(5,:)]];
            end
            % points-loop two more
            for k = nSlices-1:nSlices
                % build upo the points                
                obj.p(:,nPts*(k-1)+1:nPts*k) = ...
                    [grid2d.p;z(k)*ones(1,nPts)];  
                                       
            end
            obj.e =  [obj.e [grid2d.t(1:3,:)+nPts*(nSlices-1);...
               zeros(1,size(grid2d.t,2));...
               grid2d.t(4,:)+max(obj.e(5,:))]];
            
            
            % local block definition function, different from tetrahedra
            function block = expandFirstSlice(t,nPts)
                block = zeros(6,size(t,2));
                for kb = 1:size(t,2)
                    block(:,kb) =  [t(1,kb);...
                        t(2,kb);...
                        t(3,kb);...                        
                        nPts+[t(1,kb);...
                                t(2,kb);...
                                t(3,kb)]];
                end
            end  
            
            function block = expandBoundary(e,nPts,nBdConds)
                % expands the boundary conditions to the mantle 
                % of the cylinder
                    block  = [e(1:2,:);e(2:-1:1,:)+nPts];
                    block = [block;e(5,:)+nBdConds];
                 
            end
        end  
    end
    
    methods(Static) 
        function[qval,gval,hval,rval] = bCoefficients(p,b)
            % compute the boundary coefficients
            % 
            % q,g,h,r are srings defining everything that can be 
            % evaluated by eval. The
            % independent variables must be named by x, y, z
            % 
            % are three or four points    
            x = sum(p(1,:))/size(p,2); %#ok<*NASGU>
            y = sum(p(2,:))/size(p,2);
            z = sum(p(3,:))/size(p,2);    


            m = b(2);
            qval = 0;
            gval = 0;
            hval = 0;
            rval = 0;
            lengthq = b(3);
            lengthg = b(4);
            if m == 0 % only Neumann BCs        
                qval = eval(char(b(5:5+lengthq-1)));    
                gval = eval(char(b(5+lengthq:5+lengthq+lengthg-1)));
            else % only Dirichlet BCs
                lengthh = b(5);
                lengthr = b(6);
                char(b(9:9+lengthh-1));        
                hval = eval(char(b(9:9+lengthh-1)));
                rval = eval(char(b(9+lengthh:9+lengthh+lengthr-1)));
            end  
        end
    end      
end

