classdef grid1D < gridd
     properties(SetAccess = protected)
        nPointsInElements = 2;
    end
    
    properties(Constant = true)
        spaceDimension = 1;
    end
    methods(Access = public)       
        function obj = grid1D()
            switch nargin
                case 0
                    % empty object                                    
                otherwise
                    throw(obj.wrongNumberInputs);
            end
            obj.nPointsInElements = 2;
        end
        
        function diam = h(obj)
            % diam = obj.h() 
            % Computes the meshsize h aka (in 2D)
            % the(outer) diameter of all triangles
            diam = obj.p(2:end)-obj.p(1:end-1); 
        end
        
        function b = isBoundary(obj)
            %  b = isBoundary(gt) 
            % Method that indicates the boundary points.
            % b is a logical vector of length(#triangles)                          
            b = false(1,size(obj.p,2));        
            b(1)= true;
            b(size(obj.p,2)) = true;
        end  
        
        function refineUniformly(obj,refines)
            switch nargin
                case 1
                    n = 1;
                case 2
                    n = refines;
                otherwise
                    obj.wrongNumberInputs.throwAsCaller;
            end
            for k = 1:n
                obj.refineMesh;
            end
        end
        
        function refineMesh(obj,varargin)
            % Methods for uniform refinement of grid
            % obj.refineMesh();
            % obj.refineMesh(elements);
            % elements is an index vector of to refined elements
            if obj.isExtended
                % Brute force re-initialization and re-extension
                x = obj.p(1:obj.ngpts);
                obj2 = grid1D();
                obj2.interval(x);                 
                obj2.refineMesh(varargin{:});                
                obj2.extendMesh; 
                
                obj.p = obj2.p;
                obj.e = obj2.e;
                obj.t = obj2.t;

                obj.b = obj2.b;
                obj.ngpts = obj2.ngpts;
                obj.isExtended = obj2.isExtended;
                obj.nPointsInElements = obj2.nPointsInElements;
                obj.nPointsInElements;
                return
            end
            switch nargin
                case 1 % Uniform refinement of all elements  
                    x = obj.p(1,:);
                    x = unique([x,obj.midpts]);
                    obj.interval(x); 
                case 2
                    elements = varargin{1};
                    if ~(isnumeric(elements)&&max(elements)<=obj.nElements...
                            &&min(elements)>=1)
                        obj.wrongFormat.throw;
                    else
                        x = obj.p(1,:);                      
                        x2 = 0.5*(obj.p(1,elements(1:end))+obj.p(1,elements(1:end)+1));
                        x = unique([x,x2]);
                        obj.interval(x);
                    end
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end
               
        function mpt = midpts(obj)
             mpt = 0.5*(obj.p(1:end-1)+ obj.p(2:end));
        end
        
        function plot(obj,varargin)
            % plot method for class gridd
            % obj.plot()
            switch nargin
                case 1
                    hold on 
                    if obj.nPoints < 25
                        plot(obj.p,zeros(size(obj.p,2)),'.');                        
                    else
                        plot(obj.p,zeros(size(obj.p,2)),':');
                    end
                case 2
                    hold on
               
                    if obj.nPoints < 25
                        plot(obj.p,varargin{1},'-');                        
                    else
                        plot(obj.p,varargin{1},'-');
                    end                     
                     
                otherwise
                    hold on
                    plot(obj.p,varargin{1},varargin{2:end});
            end
        end
        
        function obj2 = copy(obj1)
            % hard copy method
            % obj2 = obj1.copy()
            obj2 = grid1D();         
            obj2.p = obj1.p;
            obj2.e = obj1.e;
            obj2.t = obj1.t;
            obj2.b = obj1.b;
            obj2.ngpts = obj1.ngpts;
            obj2.isExtended = obj1.isExtended;
            obj2.nPointsInElements = obj1.nPointsInElements;
        end
          
        
        function initMesh(obj,varargin)
              warning(['The call of the initMesh ',...
                'is obsolate and will be removed in future releases.',...
                ' Use instead interval with the same arguemnts to greate 1D grids.'])
            obj.interval(varargin{:});
        end
        
         
        
        
        function interval(obj,geo,hmax)
           
            switch nargin
                case 2 % obj and geometry
                    if ~isa(geo,'double')
                        throw(obj.wrongClass)
                    end
                    obj.p = geo(:)'; % row vector
                    obj.t = [1:size(geo,2)-1 
                             2:size(geo,2)];
                case 3                                    
                    try
                        geo = geo(:)';
                        if isa(hmax,'char')% hmax could be a char
                            hmax = eval(hmax); 
                        end                        
                        obj.p = geo; 
     nx=round((geo(2)-geo(1))/hmax/2)+1; % HU: start with mesh size 2*hmax
     obj.p=linspace(geo(1),geo(2),nx); % (1 mesh-ref somehow needed ...)
                        h = max(obj.p(2:end)-obj.p(1:end-1));
                        if h > hmax
                            while h > hmax
                                % refine
                                obj.refineMesh; break; %HU
                                h = max(obj.p(2:end)-obj.p(1:end-1));                                
                            end
                        end
                    catch ME
                        ME.throwAsCaller;
                    end                    
                otherwise
                   obj.wrongNumberInputs.throw;
            end
            % pro forma the edge field 
            obj.e = [1 obj.nPoints
                0 0
                0 0
                0 0
                1 2
                0 0];   
        end

        
        
        function bval = convCoefficientsMpt(obj,b)
            % ccoefficients static method to evaluate the b vector
           x = obj.p; 
            if isa(b,'function_handle') || isa(b,'inline')
                bval = feval(b,p);
            elseif isa(b,'char')
                bval = eval(b).*ones(1,n);
            elseif isa(b,'numeric')||isa(b,'deriv')
                if length(b)==obj.nPoints                    
                    % f vektor
                    bval = obj.point2Center(b);
                elseif length(b)==obj.nElements
                    bval = b;    
                elseif length(b)==1,
                    % skalar
                    bval = b*ones(1,obj.nElements);
                else
                    error('wrong sized b')
                end
            elseif isa(b,'inline')
                bval = b(p); 
                bval = 0.5*(bval(2:end)+bval(1:end-1));
            else
                error('wrong formated b')
            end
            % value at the midpts
            
        end
        
        function [cval,aval,fval] = aCoefficientsMpt(obj,c,a,f)
 
            n = obj.nElements;
 
            x = obj.point2Center(obj.p);
            if isa(c,'function_handle') || isa(c,'inline')
                cval = feval(c,x);
            elseif isa(c,'char'),
                cval = eval(c).*ones(1,n);
            elseif isa(c,'numeric')
                if length(c)==n
                    % f vektor
                    cval = c(:)';
                elseif length(c)==obj.nPoints
                    cval = obj.point2Center(c);
                    
                elseif length(c)==1,
                    % skalar
                    cval = c*ones(1,n);
                else
                    error('wrong sized c')
                end
            elseif isa(c,'inline')
                cval = c(p);
            else
                error('wrong formated c')
            end
            
            if isa(a,'function_handle') || isa(a,'inline')
                aval = feval(a,p);
                
            elseif isa(a,'char'),
                aval = eval(a).*ones(1,n);
                %ones(n,1);
            elseif isa(a,'numeric')
                if length(a)==n
                    % f vektor
                    aval = a; 
                    aval = aval(:)';
                elseif length(a)==obj.nPoints
                    aval = obj.point2Center(a);
                elseif length(a)==1,
                    % skalar
                    aval = a*ones(1,n);
                else
                    error('wrong sized a')
                end
            elseif isa(a,'inline')
                aval = a(p);
            else
                error('wrong formated a')
            end
            
            if isa(f,'function_handle') || isa(f,'inline')
                fval = feval(f,x);
            elseif isa(f,'char'),
                fval = eval(f).*ones(1,n);
            elseif isa(f,'numeric')
                if length(f)==n
                    % f vektor
                    fval = f;
                elseif length(f)==obj.nPoints
                    fval = obj.point2Center(f);     
                elseif length(f)==1,
                    % skalar
                    fval = f*ones(1,n);
                else
                    error('wrong sized f')
                end
            elseif isa(f,'inline')
                fval = f(p);
            else
                error('wrong formated f')
            end
        end
                       
        function element = pointToElementIndex(obj,pt)
            [~,indx] = sort(abs(obj.p-pt));
            element = indx(1);
        end
        
        function  dy = gradient(obj,f,varargin)
            dy = (f(2:end)-f(1:end-1))'...
                ./(obj.p(1,2:end)-obj.p(1,1:end-1));
        end
        
        
        function extendMesh(obj)
            % obj.extendMesh
            % Extends the mesh to use it for hig
            % Create a mesh where the additional point integratet into the
            % t structure!!! 
            % Idea: Compute the points on the midle of every subinterval
            % (Element) Add this to the Element data. An element is defined
            % by all points. Edges leave untoughed since we are in 1D
            if obj.isExtended
                return
            end
            % midpoints Note the order of operations is importand to work.
            xad = 0.5*(obj.p(2:end)+obj.p(1:end-1));               
            obj.t= [obj.t(1,:);...
                obj.nPoints+(1:length(xad));...
                obj.t(2,:)];    
            obj.ngpts = obj.nPoints;
            obj.p = [obj.p,xad];  
            obj.isExtended = true;    
            obj.nPointsInElements = 3;       
        end
        
    end  % public block ends
    
    methods(Access = private)
        function N = neighbours(obj,varargin)            
            % N = neigbours(gt[,indx]) 
            % computes to every triangle the 
            % index of neighbored triangles              
            nt = length(obj.t(1,:));            
            if nargin==2,
                indx = varargin{1};
            else
                indx = 1:nt;
            end            
            N = sparse_null(3,nt);
            
            for k = 1:length(indx),
                nk = obj.ent(indx(k));
                [i] = find(nk==k);
                nk(i)=[];
                N(1:length(nk),indx(k))=nk;
            end
        end  
    end
        
    methods(Hidden)
        % all methods from handle
        % intended not to shown...
        % addlistener(obj);  ge(obj)
        % le(obj)
        % ne(obj)
        % notify(obj)
        % gt(obj)
        % eq(obj)
        % findobj(obj)
        % findprop(obj)
        % lt(obj)        
    end
    
    methods(Access = public)           
        function [sidelength,area] = sideLengthAndArea(obj)
        %Method to compute  side lengths and areas of elements.
        % Results are vectors of length nElements
        % (c) 2015 by Uwe Prüfert
            % Let them one 
            sidelength = ones(1,obj.nElements);
            % semi-last line minus first line of t are the left and right
            % points of the sub-intervals
            %area =  obj.p(:,obj.t(end-1,:))-obj.p(:,obj.t(1,:)); 
            area =  obj.p(:,obj.t(2,:))-obj.p(:,obj.t(1,:)); 
            % Make it  compatible to extendedMeshs
            area = area(1:obj.nElements);            
        end  
    end
    
end