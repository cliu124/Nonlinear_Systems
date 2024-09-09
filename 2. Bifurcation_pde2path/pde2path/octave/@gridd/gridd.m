classdef gridd < handle
    %properties(SetAccess = protected)
    properties(Access=public)
         % p e t = points, edges triangle, as in PDE tool defined.
        p,t,e % points, edges, triangles 
        rlong=0; % HU    
        b % boundary condition matrix
        nPointsInElements % # points in element, 2 in 1D, 3 or 4 in 2D, 4 or 6 in 3D                         
        isExtended= false % if expanded then e.g. additional points added to p
        indexOfDirichletBoundarySegments=[]; 
        indexOfRobinBoundarySegments=[]; 
    end 
    properties(SetAccess=protected)        
        ngpts % number of grid points w/o extended points
    end
    
    properties(Constant)
        wrongNumberInputs = MException('GRIDD:WRONGNUMBERINPUTS',...
            'The number of arguments is wrong, check it!');
        wrongFormat = MException('GRIDD:WRONGFORMAT',...
            'Input has the wrong format, check it.');
        wrongClass = MException('GRIDD:WRONGCLASS',...
           'Wrong argument class, check it!');
        wrongNumberPoints = MException('GRIDD:WRONGPOINTS',...
            'Wrong number of points, check it!'); 
        privateProp = MException('GRIDD:PRIVAT',...
            'The property is declared as private or do not exist.'); 
        noIndexing = MException('GRIDD:INDEX',...
            'No arrays of gridd allowed.');
        wrongNumberOfBoundaryConditions = MException(...
            'GRIDD:WRONGBOUNDARYCOND',...
            'The number of boudary conditions in b and defined by geometry don''t match.')
    end 
    
    methods(Access = public)
     %    function display(o); fprintf('gridd-display\n');  endfunction
         function E = center2PointMatrix(obj)
    idx1 = reshape(obj.t(1:obj.nPointsInElements,:),...
        1,obj.nElements*obj.nPointsInElements);
    idx2 = reshape(repmat(1:obj.nElements,obj.nPointsInElements,1),...
        1,obj.nElements*obj.nPointsInElements); 

    E = sparse(...
        idx2,...
        idx1,...
        reshape(ones(obj.nPointsInElements,1)*ones(1,obj.nElements),1,obj.nPointsInElements*obj.nElements),...
        obj.nElements,obj.nPoints)';
    E = sparse(1:obj.nPoints,1:obj.nPoints,1./sum(E,2))*E;            
         end
        function E = point2CenterMatrix(obj)     
               %     y2 = M*y1 is the same as y2 = obj.point2Center(y1)
           
            idx1 = reshape(obj.t(1:obj.nPointsInElements,:),1,obj.nElements*obj.nPointsInElements);
            idx2 = reshape(repmat(1:obj.nElements,obj.nPointsInElements,1),...
                1,obj.nElements*obj.nPointsInElements);           
            E = sparse(idx2,idx1,ones(1,obj.nPointsInElements*obj.nElements),...
                obj.nElements,obj.nPoints)/obj.nPointsInElements;            
        end     
        
        function disp(obj)
            % disp of class gridd
            % obj.disp     
            % (c) 2013 Uwe Prüfert 
            if  isempty(obj.p) 
                fprintf('Empty grid object\n');
            else
                disp('Grid object with')
                fprintf(['   ',num2str(length(obj.p(1,:))),' mesh points\n']);                               
                fprintf(['   ',num2str(length(obj.e(1,:))),' edges/lateral faces \n']);                
                fprintf(['   ',num2str(length(obj.t(1,:))),' elements.\n\n']);                
            end
        end
        
        function ne =  nElements(obj)
            % obj.nElements gives back the number of elements of the
            % mesh
            % (c) 2013 Uwe Prüfert 
            ne = size(obj.t,2);
        end
        
        function np =  nPoints(obj)
            % obj.nPoints gives back the number of points in the
            % mesh
            % (c) 2013 Uwe Prüfert 
            
                np = size(obj.p,2);
            
        end
        
        function ne =  nEdges(obj)
            % obj.nEdges gives back the number of  edges/vertices
            %  
            ne = size(obj.e,2);
        end
        
        function nb =  nBoundarySegments(obj)
            % obj.nBoundarySegments gives back the number of boundary segments
            % (c) 2013 Uwe Prüfert  
            nb = length(unique(obj.e(5,:)));
        end
        
        
        function  makeBoundaryMatrix(obj,varargin)
            % Method that combines boundary 
            % vectors to a boundary matrix. 
            % b = makeBoundaryMatrix(b1,b2,...,bn)
            % number of arguments should be the 
            % same as number of boundary segments.
            % A call with only one argument will set the boundary
            % condion to all boundary segments.
            % (c) 2013 Uwe Prüfert 
        %    obj, varargin, nargin, pause  % HU 
            switch nargin
                case 1
                    throw(obj.wrongNumberOfBoundaryConditions);
                case 2
                    % one for all
                    bound = varargin{1};
                    for k = 2:obj.nBoundarySegments 
                        k1 = size(bound,2);
                        k2 = size(varargin{1},2);
                        bound(:,k1+1:k1+k2) =  varargin{1};    
                    end
                    obj.b = bound;
                    n = size(obj.b,2);  
                otherwise
                    % every segment must have an entry
                    if obj.nBoundarySegments == nargin-1
                        bound = varargin{1};
                        for k = 2:length(varargin) 
                            k1 = size(bound,2);
                            k2 = size(varargin{k},2);
                            bound(1:size(varargin{k},1),k1+1:k1+k2)=varargin{k};    
                        end
                        obj.b = bound;
                        n = size(obj.b,2);  
                    else
                        throw(obj.wrongNumberOfBoundaryConditions);
                    end
            end                      
            if obj.nBoundarySegments~=n     
                throw(obj.wrongNumberOfBoundaryConditions);
            end
        end
        
        function b = dirichletBC(obj,h,r)
            % Short call for creating Dirichlet BC
            % Syntax:
            % b = obj.dirichletBC()     homogenious Dirichlet BC
            % b = obj.dirichletBC(h,r)  general Dirichlet BC
            % (c) 2013 by Uwe Prüfert 
            switch nargin
                case 1
                    b = obj.boundaryCondition([],[],'1','0'); 
                case 2
                    b = obj.boundaryCondition([],[],'1',h);  
                case 3
                    switch class(h)
                        case 'double'
                            h = num2str(h);
                        case 'char'
                            % ok
                        otherwise
                            throw(obj.wrongClass)
                    end
                    switch class(r)
                        case 'double'
                            h = num2str(r);
                        case 'char'
                            % ok
                        otherwise
                            throw(obj.wrongClass)
                    end
                    b = obj.boundaryCondition([],[],h,r);
                otherwise
                    throw(obj.wrongNumberInputs)
            end
        end
        
        function b = neumannBC(obj,g)
            %  Short call for creating neumann BC
            %  Syntax:
            %  grid.neumannBC(g)
            b = robinBC(obj,0,g);
        end
        
        function b = robinBC(obj,q,g)
            %  Short call for creating Robin BC
            %  Syntax:
            %  b = robinBC()    homogenious Neumann BC
            %  b = robinBC(q)   homogenious Robin BC
            %  b = robinBC(q,g) general Robin BC
            % (c) 2013 by Uwe Prüfert 
            
            switch nargin
                case 1
                    b = obj.boundaryCondition();
                case 2 
                    switch class(q)
                        case 'double'
                            q = num2str(q);
                        case 'char'
                            % ok
                        otherwise
                            throw(obj.wrongClass)
                    end
                    b = obj.boundaryCondition(q);                    
                case 3
                    switch class(q)
                        case 'double'
                            q = num2str(q);
                        case 'char'
                            % ok
                        otherwise
                            throw(obj.wrongClass)
                    end
                    switch class(g)
                        case 'double'                            
                            g = num2str(g);
                        case 'char'
                            % ok
                        otherwise
                            throw(obj.wrongClass)
                    end
                    b = obj.boundaryCondition(q,g);
                otherwise
                    throw(obj.wrongNumberInputs)
            end
        end
        

       
        function b = boundaryCondition(obj,q,g,h,r)
            % Defines boundary condition matrix in PDEtool style
            % Syntax:
            % b = obj.boundaryCondition            homogenious Neumann BCs  
            % b = obj.boundaryCondition(g)         Neumann BCs
            % b = obj.boundaryCondition(q,g)       Robin BCs 
            % b = obj.boundaryCondition([],[],h,r) Dirichlet BC
            %
            % (c) 2013 by Uwe Prüfert 
            
            % handle the "minimalistic syntax" options
            switch nargin                 
                case 1
                    q = '0';
                    g = '0';
                    h = [];
                    r = [];
                case 2  
                    g =  q;
                    q = '0';
                    h = [];
                    r = [];
                case 3 
                    h = [];
                    r = [];    
                case 5 
                    if ~(((isempty(g)&&isempty(q))&& ~(isempty(h)&&isempty(r)))||...
                            (~(isempty(g)&&isempty(q))&& (isempty(h)&&isempty(r))))
                        ME = MException('obj:WRONGDEFINITION',...
                            'You can define only Dirichlet OR Robin BCs.');
                        throw(ME);
                    end
                otherwise
                    throw(obj.wrongNumberInputs)
            end
            if isa(g,'double')
                g = num2str(g);
            end
            if isa(q,'double')
                q = num2str(q);
            end
            if isa(r,'double')
                r = num2str(r);
            end
            if isa(h,'double')
                h  = num2str(h);
            end
            if isempty(h)&&isempty(r),
                b = 1;                       % only Neumann BCs
                b(2,1) = 0;
                b(3,1) = length(q);
                b(4,1) = length(g);
                b=[b;double(q)';double(g)'];
            elseif isempty(q)&&isempty(g)   % Dirichlet BCs, setting the Neumann part  
                b = 1;                       % to zero
                b(2,1) = 1;
                b(3,1) = 1;
                b(4,1) = 1;
                b(5,1) = length(h);
                b(6,1) = length(r);
                b(7,1) = 48;
                b(8,1) = 48;
                b = [b;double(h)';double(r)'];    
            end
        end
        
        
        function val = point2Center(obj,y)
            % point2Center - ecaluates y given at the nodes of a
            % triangulation at the center of elements
            % (c) 2014 Uwe Prüfert
            y = y(:); 
 
            switch obj.nPointsInElements
                case 2                       
                   val = sum([y(obj.t(1,:))' 
                       y(obj.t(2,:))'])/2;
                case 3                    
                   val = sum([y(obj.t(1,:))'
                            y(obj.t(2,:))'
                            y(obj.t(3,:))'])/3;
                case 4                    
                   val = sum([y(obj.t(1,:))'
                            y(obj.t(2,:))'
                            y(obj.t(3,:))'
                            y(obj.t(4,:))'])/4;
                case 6
                   
                   val = sum([y(obj.t(1,:))'
                            y(obj.t(2,:))'
                            y(obj.t(3,:))'
                            y(obj.t(4,:))'
                            y(obj.t(5,:))'
                            y(obj.t(6,:))'])/6;
                otherwise
                    MException('PDE:point2Center:UnexpectedNumber',...
                        'Unexpected number of points in element').throw;
            end
        end 
        
        function varargout = nearestPointInGrid(obj,points)
            %%% 
            
            indx = zeros(size(points,2),1);
            val = zeros(size(points,2),1);
            for k = 1:size(points,2)
                [val(k),indx(k)] = min(sum((obj.p-points(:,k)*ones(1,obj.nPoints)).^2));                 
            end
            switch nargout                
                case 1
                    varargout{1} = indx;
                case 2
                    varargout{1} = val;
                    varargout{2} = indx;
                otherwise
                    % ignore
            end
        end
        
        
    end
    methods(Access = public, Hidden)
        function setPET(obj,p,e,t)
            obj.p = p;
            obj.e = e;
            obj.t = t;
        end
    end
end