classdef pde < handle 
    properties(Constant,Access=protected)
        % some error objects
        wrongClassID = 'PDE:wrongClass';
        wrongClassStr = 'Wrong argument class, it must be ';
        wrongNumberArgumentsID =  'PDE:wrongNumberArguments';
        wrongNumberArgumentsStr = 'Wrong number of arguments, it must be ';
        wrongSizedArgumentID = 'PDE:wrongSizeArguments';
        wrongSizedArgumentStr = 'Wrong size of Arguments, it must be '; 
        invalidOperationID = 'PDE:invalidOperation';
        invalidOperationStr = 'The operation is not valid: ';
        unknownPropertiesID = 'PDE:unknownField';
        unknownPropertiesStr = 'Unknown property:';        
        wrongInputFormatID = 'PDE:wrongInputFormat';
        wrongInputFormatStr = 'Wrong input format.';
        
        notInitialized = MException('PDE:notInitialized',...
            'Object should be initialized.');
        notAllowed = MException('PDE:notAllowed',...
            'This operation is not allowed for class pde.'); 
    end
    properties(SetAccess = public)
        fem      % fem-object        
        grid     % grid-object        
        time     % [tt_start t_end]
        y;       % solution vector or matrix        
    end 
    
    properties(SetAccess = public, Hidden)
        odeOptions= odeset('Stats','off');
        solverOptions= struct('solver','gauss','solverTol',1e-6);
    end
    
    properties(Access = public)        
        pattern % sparsity pattern for Jacobian  
        initialized= false;  
        % All the FEM matrices for standard PDEs
        K,M,C,H % Stiffness,  Mass, Convection, Dirchechlet BC      
        Q,F,G,R % Neumann BC, RHS, RHS Neumann, RHS Dirichlet
        A,b,D % general linear matrix, general rhs, Integral matrix for d*dy/dt       
    end
    
    properties(Access = protected,Hidden)    
        g,fc  % helper props for numjac
    end
    
    methods(Static)
	function indx = selectElements2Refine(E,sigma)
            % Element selection based on Dörflers method
            % Should be overwritten for systems
            % (c) 2015 Uwe Prüfert 
            switch nargin
                case 1
                    sigma = 0.5;
                case 2
                    %
                otherwise
                    MException(obj.wrongNumberArgumentsID,...
                         [obj.wrongNumberArgumentsStr,' 1 or 2.'])
            end
            [vals,idx] = sort(E,'descend');
            sumeta=cumsum(vals);
            ntref = find(sumeta>=sumeta(end)*sigma,1);
            indx=idx(1:ntref); 
    end  
    end
    
    
    methods(Access = public) 
        function setBoundaryConditions(obj,varargin)
        %% setBoundaryConditions
            % 
            % Sets the boundary conditions.
            %   
            %   pde.setBoundaryConditions(cond1,...,condN)
            % 
            % cond1 .. condN must be pairs type,value
            % 
            % Type can be 'Dirichlet', 'Neumann' or 'Robin'
            % Value can be char, or double. 
            %
            % Value must be for types 'Dirichlet' and 'Neumann' a scalar char or
            % double and for type 'Robin' a vector of length two of double or 
            % a cell array of length two for char.
            %
            % If there is only one boundary conditon on all boudry segments
            % it is sufficient to define one boundary condition.
            %
            % Examples
            %
            % Homogenious Neuman on all boundary segments.
            %
            %     pde.setBoundaryConditions(...
            %         'Neumann','0');
            % 
            % Two boundary segments. No1 homogenious Dirichlet, no2
            % homogenious Neumann boundary condition.
            % 
            %     pde.setBoundaryConditions(...        
            %         'Dirichlet','0',...
            %         'Neumann','0');
            % 
            % One every boundary segment Dirichlet $u = cos(s)$, where $s$ is
            % the nomalized arclenght parameter (2D domains only).
            % 
            %     pde.setBoundaryConditions(...        
            %         'Dirichlet','cos(s)');
            %
            % 
            % Robin type
            % 
            % $$\vec n \cdot  c \nabla u + u = cos(s)$$
            %
            %     pde.setBoundaryConditions(...        
            %         'Robin',{'1','cos(s)'});
            
            bcs = cell(1,ceil(length(obj.grid)/2));
            l = 1;
            for k = 1:2:length(varargin)                
                switch varargin{k}
                    case 'Dirichlet'
                        if length(varargin{k+1})~=1
                           MException('pde:setBoundaryCondition:missformedArgument' ,...
                                'The argument  must be scalar').throwAsCaller;    
                        end
                        bcs{l} = obj.grid.dirichletBC(varargin{k+1});
                        l = l+1;
                    case 'Neumann'
                        if length(varargin{k+1})~=1
                           MException('pde:setBoundaryCondition:missformedArgument' ,...
                                'The argument  must be scalar').throwAsCaller;    
                        end
                        bcs{l} = obj.grid.neumannBC(varargin{k+1});
                        l = l+1;
                    case 'Robin'
                        if length(varargin{k+1})~=2
                           MException('pde:setBoundaryCondition:missformedArgument' ,...
                                'The argument  must be a double or cell vector of lenght two.').throwAsCaller;    
                        end
                        switch class(varargin{k+1})
                            case 'cell'
                                bcs{l} = obj.grid.robinBC(varargin{k+1}{:});
                            case 'double'
                                bcs{l} = obj.grid.robinBC(varargin{k+1}(1),varargin{k+1}(1));
                            otherwise
                        end
                        l = l+1;
                    otherwise
                        MException('pde:setBoundaryCondition:unknownArgument',...
                            ['The argument ''',varargin{k},''' is not known.',...
                            ' Allowed values are ''Dirichlet'',',...
                            ' ''Neumann'' or '' Robin''']).throwAsCaller;
                end
            end            
            obj.grid.makeBoundaryMatrix(bcs{:});
    end
        
        function [varargout] = gradient(obj,u,mode)
    %% gradient
    %
    % * gradient  IN u:double[,mode:char]
    %            OUT dx:double[, dy:double[ dz:double]]  
    %
    % Method that computes the gradient of u
    %
    %    [dx[,dy[,dz]]] = gradient(grid,u[,mode])
    %
    % u must be a vector of length number of grid points. mode can
    % be 'points' or 'elements'. It controls if the gradient is
    % evaluated at the center of the elements or at the grid
    % points. If no mode is given, 'elements' will be used.
    % 
    % It gives back row vectors of length nElement/nPoints
    % with the values of            
    % $\nabla u$
    % . The number of output arguments depends on the dimension of
    % the spatial domain. 
    % Examples
    % 
    % *1D domains*
    %
    %    dx = pde.gradient(u)
    %
    % *2D domains*
    %
    %    [dx,dy] = pde.gradient(u)
    %
    % If only one derivative is of interest, call e.g. for dy
    %
    %    [~,dy] = pde.gradient(u);
    % 
    % For dx call it with only one argument. The calls
    % 
    %    [dx,~] = pde.gradient(u);
    %    dx = pde.gradient(u);
    %
    % are equivalent.
    %
    % *3D domains*
    %
    %    [dx,dy,dy] = pde.gradient(u)
    %
    % If only one derivative is of interest, call e.g. for dy
    %
    %    [~,dy,~] = pde.gradient(u);
    %
    % or
    %
    %    [~,dy] = pde.gradient(u);
    %
    % For dx the calls
    % 
    %    [dx,~,~] = pde.gradient(u);
    %    [dx,~] = pde.gradient(u);
    %    dx = pde.gradient(u);
    %
    % are equivalent.



    if nargin ==2
        mode = 'elements';
    end
    if length(u) == obj.grid.nPoints  
        switch obj.grid.spaceDimension                   
            case 1
                DX = obj.fem.gradientMatrices(obj.grid);
                switch mode
                    case 'points'                              
                        varargout{1} = (obj.grid.center2PointMatrix*(DX*u(:)))';
                    case 'elements'
                        varargout{1} = (DX*u(:))';
                    otherwise
                        MException('pde:gradient:wrongmode',...
                            'The mode must be *points* or *elements*').throwAsCaller;
                end
            case 2
                [DX,DY] = obj.fem.gradientMatrices(obj.grid);                         
                switch mode
                    case 'points'                               
                        varargout{1} = (obj.grid.center2PointMatrix*(DX*u(:)))';
                        varargout{2} = (obj.grid.center2PointMatrix*(DY*u(:)))';
                    case 'elements'                                 
                        varargout{1} = (DX*u(:))';
                        varargout{2} = (DY*u(:))'; 
                otherwise
                        % wrong mode error
                end
            case 3
                [DX,DY,DZ] = obj.fem.gradientMatrices(obj.grid);
                switch mode
                    case 'points'                               
                        varargout{1} = (obj.grid.center2PointMatrix*(DX*u(:)))';
                        varargout{2} = (obj.grid.center2PointMatrix*(DY*u(:)))'; 
                        varargout{3} = (obj.grid.center2PointMatrix*(DZ*u(:)))'; 
                    case 'elements'
                        varargout{1} = (DX*u(:))';
                        varargout{2} = (DY*u(:))'; 
                        varargout{3} = (DZ*u(:))'; 
                otherwise
                        % wrong mode error
                end
            otherwise
        end                
    else 
        % wrong format error
        MException('pde:gradient:wrongformat',...
            'Input must be a vector of length *no of points*.').throwAsCaller;
    end            
end
       
        function E = errorInd(obj,y,varargin)
            E = obj.fem.errorInd(obj.grid,y,varargin{:}); 
        end
        
        function initialize(obj,varargin)            
            switch nargin
                case 5  
                    [obj.K,obj.M,obj.F] = obj.fem.assema(obj.grid,...
                        varargin{1},varargin{3},varargin{4});
                    [obj.Q,obj.G,obj.H,obj.R] = obj.fem.assemb(obj.grid);
                    obj.C = obj.fem.convection(obj.grid,varargin{2});
                    obj.initialized = true;
                case 6
                    [obj.K,obj.M,obj.F] = obj.fem.assema(obj.grid,...
                        varargin{2},varargin{4},varargin{5});
                    [obj.Q,obj.G,obj.H,obj.R] = obj.fem.assemb(obj.grid);                 
                    obj.C = obj.fem.convection(obj.grid,varargin{3});
                    [~,obj.D,~] = obj.fem.assema(obj.grid,0,varargin{1},0);
                    obj.initialized = true;
                otherwise
                    MException(obj.wrongNumberArgumentsID,...
                     [obj.wrongNumberArgumentsStr,...
                     ' 5 (stationary case) or 6 (transient case)']).throwAsCaller;
            end            
        end
        
        function solve(obj,solver)
            % universal solver method
            % obj.solve(solver)
            % where solver can be on of the
            % strings:
            % EULERI, EULERILIN, BDF2, ODE15S, ODE23S, LINEAR, INEARGAUSS
            if ~ischar(solver)
                MException(obj.wrongClassID,...
                    [obj.wrongClassStr,' a string']).throwAsCaller;
            end
            switch solver
                case 'EULERI'
                    obj.euleri;
                case 'EULERIT'
                    obj.eulerit;
                case 'BDF2'
                    obj.bdf2;
                case 'BDF2STEPSIZE'
                    obj.bdf2sc;  
                case 'ODE15S'
                    obj.ode15s;
                case 'ODE23S'
                    obj.ode23s;
                case 'EULERILIN'
                    obj.eulerilin;
                case 'STEADYSTATE'
                    obj.steadystate;
                case 'STATIONARY'
                    obj.stationary
                case 'LINEAR'
                    obj.linear
                case 'LINEARGAUSS'
                    obj.linearGauss
                otherwise
                     MException(obj.invalidOperationID,...
                         [obj.invalidOperationStr,' solver',solver,...
                         ' is not supported.']).throwAsCaller; 
                    
            end
        end
        
        function K = stiff(obj)
            if ~isempty(obj.grid)                
                [K,~,~] = obj.fem.assema(obj.grid,1,0,0); 
            else
                obj.notInitialized.throwAsCaller;
            end            
        end
        
        function M = mass(obj,varargin)
            if ~isempty(obj.grid)                
                [~,M,~] = obj.fem.assema(obj.grid,0,1,0); 
            else
                obj.notInitialized.throwAsCaller;
            end
            switch length(varargin)
                case 0
                    % do nothing
                case 1                     
                    if strcmp(varargin{1},'lumped')
                        M = diag(sum(M));
                    else
                        MException(obj.wrongClassID,...
                            [obj.wrongClassStr,...
                            ' string. Only allowed value is '' lumped''']).throwAsCaller;
                    end
                otherwise
                    MException(obj.wrongNumberArgumentsID,...
                            [obj.wrongNumberArgumentsStr,...
                            ' zero or one']).throwAsCaller;
            end
        end
        
        function val = L2(obj,arg)
            % We need an initialized object
            if ~obj.initialized
                obj.notInitialized.throwAsCaller;
            end
            switch nargin
                case 1
                    % compute normu of obj.y   
                    if isempty(obj.y)
                        MException(obj.invalidOperationID,...
                            [obj.invalidOperationStr,...
                            ' L2-Norm needs a non-empty obj.y.']);
                    end
                    
                    if isempty(obj.time)
                        % no time integration
                        val = sqrt(obj.y(:)'*obj.mass*obj.y(:));
                    else
                        % integrate in time
                        val = 0;                         
                        dt = [obj.time(2)-obj.time(1),...
                            0.5*(obj.time(3:end)-obj.time(1:end-2)),...
                            obj.time(end)-obj.time(end-1)];
                        for k = 1:length(obj.time)
                            val = val+(dt(k)*sqrt(obj.y(:,k)'*obj.mass*obj.y(:,k)));
                        end
                    end
                case 2
                    if isvector(arg)&& obj.grid.nPoints~=length(arg)
                         MException(obj.wrongInputFormatID,...
                             obj.wrongInputFormatStr).throwAsCaller;
                    else
                        %
                        [xdim,tdim] = size(arg);
                        if xdim~=obj.grid.nPoints&&...
                                tdim~=length(obj.time)
                            MException(obj.wrongSizedArgumentID,...
                                [obj.wrongSizedArgumentStr,...
                                'Numbers of columns = nPoints and numbers'...
                                ' of rows = number of time steps'])
                        end
                    end
                        
                    % use arg to compute norm
                    if isempty(obj.time)||isvector(arg)
                        % no time integration
                        % arg mus be a vector
                        val = sqrt(arg(:)'*obj.mass*arg(:));
                    else
                        % integrate in time
                        val = 0;  
                        
                        
                        
                        dt = [obj.time(2)-obj.time(1),...
                            0.5*(obj.time(3:end)-obj.time(1:end-2)),...
                            obj.time(end)-obj.time(end-1)];
                        for k = 1:length(obj.time)
                            val = val+(dt(k)*sqrt(arg(:,k)'*obj.mass*arg(:,k)));
                        end
                    end
                otherwise
                    MException(obj. wrongNumberArgumentsID,...
                    [obj.wrongNumberArgumentsStr,...
                    ' zero or one']).throwAsCaller;
            end
        end 
        
        function T = evaluateAtPoint(obj,point)
            % evaluateAtPoint evaluat obj.y at single point [x,y]

            indx = obj.grid.pointToElementIndex(point);
            indx2 = unique(obj.grid.t(1:end,indx));
            T = sum(obj.y(indx2,:))./length(indx2);
        end
        
        function refineMesh(obj,indx)
            % refineMesh refines mesh and interpolates current solution to
            % the new mesh.
            % Optional argument is a vector containing the elements to
            % refine.
            % Note that in the case of transient problems, all solutions
            % were transformed to the new mesh.
            % Be aware, that refineMesh is intentionally not for updating
            % the linear system and will do nothing to force an update of
            % it.
            % (c) 2015 Uwe Prüfert
            
            
            if isempty(obj.y)
                MException('PDE:NOTVALID',...
                    ['No data in obj.y. pde.refineMesh is for ',...
                    ' refinement of mesh + data.\n',...
                    'Use pde.grid.refineMesh instead']).throw;
            end
            
            % The "essentials" of the old mesh
            nPointsOld = obj.grid.nPoints;
            p = obj.grid.p;
            
            % dimension of the domain
            dim = size(obj.grid.p,1);
            % number of elements in y / nPoints
            dimSystem = size(obj.y,1)/obj.grid.nPoints;
            switch nargin
                case 1 % w/o index                    
                    obj.grid.refineMesh()                     
                case 2
                    if dimSystem==3
                        MException('PDE:NOTIMPLEMENTED',...
                            'Local refinement is still not implemented').throw;
                    end
                    obj.grid.refineMesh(indx)                    
                otherwise                    
                    % wrong call, ME later
            end 
            % interpolate solutions           
            % Note: grid.nPoints is now the new number of points! 
            
            ynew = zeros(dimSystem*obj.grid.nPoints,size(obj.y,2));
            for k = 0:dimSystem-1  
                switch dim
                    case 1
                        ynew(k*obj.grid.nPoints+1:(k+1)*obj.grid.nPoints,:) =...
                            interp1(p,obj.y(k*nPointsOld+1:(k+1)*nPointsOld,:),...
                            obj.grid.p,'linear');
                    case {2,3}
                        % need time loop. It 
                        for l=1:size(obj.y,2)
                            interpolant = scatteredInterpolant(...
                                p',obj.y(k*nPointsOld+1:(k+1)*nPointsOld,l));
                            ynew(k*obj.grid.nPoints+1:(k+1)*obj.grid.nPoints,l) = ...
                                interpolant(obj.grid.p');
                        end
                    case 3
                         
                    otherwise
                        % ME werfen...
                end
            end
            obj.y = ynew;
            
            % We must remove the old numjac f and g!
            obj.fc = [];
            obj.g = [];
            obj.initialized = false;
        end
          
    end
    
    methods(Hidden,Access = public)
        function obj = subsasgn(obj,index,val)
            % overwrite the builtin to handle an input check...
            switch index(1).type
                case '.'
                    switch index(1).subs
                        case 'time'
                            if isvector(val)
                                obj.time = val(:)';
                            else
                                MException(obj.wrongInputFormatID,...
                                    [obj.wrongInputFormatStr,...
                                    ' Must be a double vector']).throwAsCaller;
                            end                                                   
                        case 'y'
                            if isvector(val)
                                obj.y = val(:);                                
                            else % double matrix
                                [m,n] = size(val);
                                if m==obj.grid.nPoints && n==length(obj.time)
                                    obj.y = val;
                                else 
                                    MException(obj.wrongInputFormatID,...
                                        [obj.wrongInputFormatStr,...
                                        ' Data may be wrong formated.']).throwAsCaller;
                                end
                            end
                        otherwise
                            % all common calls give through to the built-in
                            try 
                                obj = builtin('subsasgn',obj,index,val);
                            catch ME
                                ME.throwAsCaller;
                            end
                    end
                otherwise
                    % all common calls give through to the built-in
                    try 
                        obj = builtin('subsasgn',obj,index,val);
                    catch ME
                        ME.throwAsCaller;
                    end
            end
        end
    end
    
    methods(Access = protected)
        function J = jacobian(obj,t,y)
            %jacobian Computes the Jacobian by finite differences.
            % Overwrite it if an exact Jacobian is available.
            % If a sparsity pattern is given, jacobian will use it.
            % Signature:
            %    j = obj.jacobian(t,y)
            % Note: When overwritten in child classes,
            % it must have Access = protected.
            % (c) 2013 by Uwe Prüfert
           tr = 1e-6;
           thresh = tr(ones(size(y)));           
            if isempty(obj.pattern)
                [J,obj.fc] = numjac(@obj.df,t,y,obj.df(t,y),thresh,obj.fc,0);
                % Make it sparse
                J = sparse(J);
            else
                [J,obj.fc,obj.g] = numjac(@obj.df,t,y,obj.df(t,y),...
                    thresh,obj.fc,0,obj.pattern,obj.g);
            end
        end
        
        function [value,isterminal,direction] = eventfun(obj,t,y)
            %eventfun event function that stopps the evaluation of the time
            % integration. Criteria for being stationary is 
            % || dY/dt(t,Y)|| < 1e-3.
            % Overwrite if you want to have a different criteria to stepp
            % the integration.
            % Signature:
            %     [value,isterminal,direction] = eventfun(obj,t,y)
            % When overwritten in child classes,
            % it must have Access = protected.
            % See also help odeset for more informations.
            % (c) 2013 by Uwe Prüfert         
             
            value =  norm(obj.df(t,y),Inf)-obj.solverOptions.solverTol;%  
            isterminal = 1;   % stopps the integration
            direction = 0;   % both directions
        end
        
        function ode15s(obj)
            %ODE15S solves the time dependent PDE problem 
            % d*du/dt = df(t,x,u(t,x)) by using MATLAB ode15s
            % function.
            % If no matrix D is found, a Mass with d = 1 will be used.
            % obj.df must be defined properly. 
            % obj.time must be a non-empty double vector of length >2.
            % It will be overwritten by the ode solver.
            % ode15  with standard options is used for time integration.
            % If obj.pattern is non empty, a sparsity pattern for the
            % jacobian is used. The jacobian will be computed by a call to
            % obj.jacobian. The standard jacobiaon is numjac. Overwrite
            % jacobian if a better one is on hand.
            % protected method. Use it  within our own methods or use
            % obj.solve('ODE15S')
            %  
            % (C) 2013 by Uwe Prüfert
            if isempty(obj.D)
                obj.D = obj.mass;
                warning('PDE:NOTIMEMASS',...
                    'No mass for time derivative given. Use Md with d = 1.')
            end
            if isempty(obj.y)||~obj.initialized
                fprintf('Not initialized or no initial value.\n')
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end  
            
            if size(obj.D,1)~=size(obj.y,1)
                MException('PDE:WRONGD',...
                    'Dimension of D do not fit the dimension of K, M, y etc.').throw;
            end
            obj.odeOptions = odeset(...
                obj.odeOptions,...
                    'Mass',obj.D,... 
                    'Jacobian',@obj.jacobian);
            if ~isempty(obj.pattern)
                obj.odeOptions = odeset(obj.odeOptions,...                    
                    'JPattern',obj.pattern);                                
            end
             
            [obj.time,obj.y] = ode15s(@obj.df,obj.time,obj.y,obj.odeOptions);
            obj.y = obj.y';
            obj.time = obj.time(:)';
             
        end
        
        function ode23s(obj)
            %ODE15S solves the time dependent PDE problem 
            % d*du/dt = df(t,x,u(t,x)) by using MATLAB ode15s
            % function.
            % If no matrix D is found, a Mass with d = 1 will be used.
            % obj.df must be defined properly. 
            % obj.time must be a non-empty double vector of length >2.
            % It will be overwritten by the ode solver.
            % ode15  with standard options is used for time integration.
            % If obj.pattern is non empty, a sparsity pattern for the
            % jacobian is used. The jacobian will be computed by a call to
            % obj.jacobian. The standard jacobiaon is numjac. Overwrite
            % jacobian if a better one is on hand.
            % protected method. Use it  within our own methods or use
            % obj.solve('ODE15S')
            %  
            % (C) 2013 by Uwe Prüfert
            if isempty(obj.D)
                obj.D = obj.mass;
                warning('PDE:NOTIMEMASS',...
                    'No mass for time derivative given. Use Md with d = 1.')
            end
            if isempty(obj.y)||~obj.initialized
                fprintf('Not initialized or no initial value.\n')
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end  
            
            if size(obj.D,1)~=size(obj.y,1)
                MException('PDE:WRONGD',...
                    'Dimension of D do not fit the dimension of K, M, y etc.').throw;
            end
            obj.odeOptions = odeset(...
                obj.odeOptions,...
                    'Mass',obj.D,... 
                    'Jacobian',@obj.jacobian);
            if ~isempty(obj.pattern)
                obj.odeOptions = odeset(obj.odeOptions,...                    
                    'JPattern',obj.pattern);                                
            end
             
            [obj.time,obj.y] = ode23s(@obj.df,obj.time,obj.y,obj.odeOptions);
            obj.y = obj.y';
            obj.time = obj.time(:)';
             
        end
        
        function euleri(obj)
            %EULERI
            % Integrates PDE in time by implizit euler scheme.             
            % An fixed time step is used, given by obj.time. If
            % obj.time contains only t0 te, stepsize is (te-t0)/100 
            % Usuage: protected method. Use it within our own methods
            % or call
            % obj.solve('EULERI')
            % 
            % (C) 2013 by Uwe Prüfert
            tic
            % "local global" focus for jacobian and  step counters
            MaxNewtonSteps  = 25;
            n_jac = 1;
             
            if isempty(obj.y)||~obj.initialized
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end 
            if isempty(obj.D)
                obj.D = obj.mass;
                warning('PDE:NOTIMEMASS',...
                    'No mass for time derivative given. Use M(d) with d = 1.')
            end
            switch length(obj.time)
                case 0
                     MException(obj.wrongInputFormatID,...
                         [obj.wrongInputFormatStr,...
                         ' obj.time is empty.']).throwAsCaller
                case 1                     
                    obj.time = linspace(0,obj.time,100);
                case 2
                    obj.time = linspace(obj.time(1),obj.time(2),100);
                otherwise
                    % obj.time is assumed as the  vector of timesteps
            end           
            dt = obj.time(2:end)-obj.time(1:end-1);
            
            % to speed up, we compute the inverse of Jacobian.
            % Inside the loop we use this matrix and it  we be updated if
            % Newton method converges "slowly".
             
            iJ = (obj.D-dt(1)*obj.jacobian(obj.time(1),obj.y(:,1)))\speye(size(obj.D));
            
            for k = 2:length(obj.time)
                % time loop, call Newton to solve the non linear problem
                % Note that the inverse ob jacobian will possibly be
                % updated inside newton, so it mus be a give through
                % parameter.
                try
                    [iJ,obj.y(:,k)] = timestep(obj,k-1,iJ,dt(k-1));
                catch ME
                    warning('ode:timestepFailed',['An error accured. Message was: ',...
                        ME.message]);
                    break
                end
            end   
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf(['PDE:EULERI:\nElapsed time ',...
                    num2str(toc),...
                    ' sec using '...
                    num2str(n_jac),' Jacobians.\n']);   
            end
            % Definition of LOCAL Function within the method
            % 
            function [iJ,y] = timestep(obj,k,iJ,dt)  
                % performs one time step
                StepFailed = false;
                
                y = obj.y(:,k); 
                t = obj.time(k+1);
                dy =  -dt*obj.df(t,y);
                
                for k1 = 1:MaxNewtonSteps                  
                    d = -iJ*dy;  
                    [~,id] = lastwarn;
                    if strcmp(id,'MATLAB:nearlySingularMatrix')|| ...
                            strcmp(id,'MATLAB:singularMatrix')
                        StepFailed = true;               
                    end
                    
                    if norm(d,Inf)  < 1e-6
                        StepFailed = false;                        
                        break
                    end                
                    y = y + d;                    
                    dy = obj.D*(y-obj.y(:,k))-dt*obj.df(t,y); 
                    
                    %  ONLY in case of slow
                    %  convergence we compute a new Jacobian,
                    %  we compute at 10-th 15-th and
                    %  20-th iteration a new one.
                    switch k1 
                        case {10 15 20} 
                            iJ = (obj.D-dt*obj.jacobian(t,y))...
                                \speye(size(obj.D));                             
                            n_jac = n_jac+1;
                        otherwise
                            % do nothing
                    end                
                end
                
                if StepFailed                     
                    MException('newton:stepFailed',...                
                        ['Unable to find a solution: ',...
                        lastwarn]).throwAsCaller;                     
                end
                if k1 == MaxNewtonSteps
                    warning('PDE:NewtonSteps',...
                        'Maximal number of Newton Steps reached.')
                end 
                
            end
        end        
        
        function eulerit(obj)
            %EULERIT
            % Integrates PDE in time by implizit euler scheme.             
            % An fixed time step is used, given by obj.time. If
            % obj.time contains only t0 te, stepsize is (te-t0)/100 
            % The linear systems will be solved by pcg. Hence we need
            % a symmetric Jacobian.
            % Usuage: protected method. Use it in your own methods or 
            % call it by                     
            % >>obj.solve('EULERIT')
            % (C) 2013 by Uwe Prüfert
            tic
            % "local global" focus for jacobian and Step counters
            MaxNewtonSteps  = 25;
            n_jac = 1;
            tol = 1e-8;
            maxIt = 200;
             
            if isempty(obj.y)||~obj.initialized
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end 
            if isempty(obj.D)                
                obj.D = obj.mass;
                warning('PDE:NOTIMEMASS',...
                    'No mass for time derivative given. Use M(d) with d = 1.')
            end
            switch length(obj.time)
                case 0
                     MException(obj.wrongInputFormatID,...
                         [obj.wrongInputFormatStr,...
                         ' obj.time is empty.']).throwAsCaller
                case 1                     
                    obj.time = linspace(0,obj.time,100);
                case 2
                    obj.time = linspace(obj.time(1),obj.time(2),100);
                otherwise
                    % obj.time is assumed as the  vector of timesteps
            end           
            dt = obj.time(2:end)-obj.time(1:end-1);
            % to speed up, we compute the inverse of Jacobian.
            % Inside the loop we use this matrix and it  we be updated if
            % Newton method converges "slowly".
            
            

            J = (obj.D-dt(1)*obj.jacobian(obj.time(1),obj.y(:,1)));
            for k = 2:length(obj.time)
                % time loop, call Newton to solve the non linear problem
                % Note that the inverse ob jacobian will possibly be
                % updated inside newton, so it mus be a give through
                % parameter.
                [J,obj.y(:,k)] = timestep(obj,k-1,J,dt(k-1));
            end   
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf(['PDE:EULERI:\nElapsed time ',...
                    num2str(toc),...
                    ' sec using '...
                    num2str(n_jac),' Jacobians.\n']);   
            end
            % Definition of LOCAL Function within the method
            % 
            function [J,y] = timestep(obj,k,J,dt) 
                % timestep perform one timestep
                StepFailed = false;
                
                y = obj.y(:,k); 
                t = obj.time(k+1);
                dy =  -dt*obj.df(t,y);
                
                for k1 = 1:MaxNewtonSteps 
                    switch obj.solverOptions.solver
                        case 'gauss'                              
                             d = J\dy;                            
                             flag = 5;
                        case 'pcg'
                            try                                
                                L = ichol(J);                               
                                [d,flag] = pcg(J,dy,tol,maxIt,L,L');
                                
                                if flag>0
                                    MException('PDE:SOLVE:SOLVERFAILED',...
                                        'Iterative solver failed').throw;
                                end
                            catch ME
                                ME.message
                                d = J\dy;                                 
                            end
                        case 'bicg'
                            try                                
                                [L,U] = ilu(J,struct('type','nofill'));
                                [d,flag] = bicg(J,dy,tol,maxIt,L,U);                               
                                if flag>0
                                    MException('PDE:SOLVE:SOLVERFAILED',...
                                        'Iterative solver failed').throw;
                                end
                            catch ME
                                ME.message
                                d = J\dy;
                                
                            end
                        otherwise
                            MException('PDE:SOLVEUNKNOWNSOLVER',...
                                ['Solver' obj.solverOptions.solver ,'is not valid.'])
                    end
                    
                    
                    [~,id] = lastwarn;
                    if strcmp(id,'MATLAB:nearlySingularMatrix')|| ...
                            strcmp(id,'MATLAB:singularMatrix')
                        StepFailed = true;               
                    end
                    
                    if norm(d,Inf)  < 1e-6
                        StepFailed = false;                        
                        break
                    end                
                    y = y - d;                    
                    dy = obj.D*(y-obj.y(:,k))-dt*obj.df(t,y); 
                    
                    %  ONLY in case of slow
                    %  convergence we compute a new Jacobian,
                    %  we compute at 10th 15 and
                    %  20th iteration a new one.
                    switch k1 
                        case {10 15 20} 
                            J = (obj.D-dt*obj.jacobian(t,y));                             
                            n_jac = n_jac+1;
                        otherwise
                            % do nothing
                    end                     
                    switch flag 
                        case 0
                            % everything is fine;-)
                        case 1  
                            fprintf('solver iterated MAXIT times but did not converge.\n');
                        case 2  
                            fprintf('preconditioner M was ill-conditioned.\n');
                        case 3  
                            fprintf('solver stagnated (two consecutive iterates were the same.\n');
                        case 4  
                            fprintf(['One of the scalar quantities calculated ',...
                                        'during iteration became too',...
                                        'small or too large to continue computing.\n']); 
                    end
                end
                
                if StepFailed                     
                    MException('newton:stepFailed',...                
                        ['Unable to find a solution: ',...
                        lastwarn]).throwAsCaller;                     
                end
                if k1 == MaxNewtonSteps
                    warning('PDE:NewtonSteps',...
                        'Maximal number of Newton Steps reached.')
                end 
                
            end
        end        
        
        function bdf2(obj)
            %BDF2
            % Integrates PDE in time by (implizit) BDF-2 scheme.             
            % A fixed time step is used, given by obj.time. If
            % obj.time contains only te, t0 = 0 is assummed.
            % If obj.time contains t0 and te, stepsize is (te-t0)/100 
            % Usuage:
            % >>obj.bdf2()
            % (C) 2013 by Uwe Prüfert
            tic
            % "local global" focus for jacobian and Step counters
            MaxNewtonSteps  = 25;
            n_jac = 1;            
            if isempty(obj.y)||~obj.initialized
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end 
            if isempty(obj.D)
                obj.D = obj.mass;
                warning('PDE:NOTIMEMASS',...
                    'No mass for time derivative given. Use M(d) with d = 1.')
            end
            switch length(obj.time)
                case 0
                     MException(obj.wrongInputFormatID,...
                         [obj.wrongInputFormatStr,...
                         ' obj.time is empty.']).throwAsCaller
                case 1                     
                    obj.time = linspace(0,obj.time,100);
                case 2
                    obj.time = linspace(obj.time(1),obj.time(2),100);
                otherwise
                    % obj.time is assumed as the  vector of timesteps
            end 
            % First, we will start with a few implicite
            % EULER steps
            % We apply solveEULERI on the time interval [t(0),t(1)] 
            
            % We must save the original time interval
            timesteps = obj.time;
            obj.time = linspace(obj.time(1),obj.time(2),5);
            
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf('PDE:BDF2: start with EULER Step\n')
            end
            obj.euleri;
            
            obj.y = obj.y(:,[1,end]); 
            obj.time = timesteps;
            dtime = timesteps(2:end)-timesteps(1:end-1);
    
            
            % to speed up, we compute the inverse of Jacobian.
            % Inside the loop we use this matrix and we be updated if
            % Newton method converges "slowly".
            
            djInv = ((2*dtime(2)+dtime(1))/(dtime(2)^2+dtime(2)*dtime(1))*obj.D...
                -obj.jacobian(obj.time(2),obj.y(:,2)))\...
                speye(size(obj.D));
            
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf('PDE:BDF2: Start BDF loop.\n')
            end
            for k = 3:length(obj.time)
                % time loop, call Newton to solve the non linear problem
                % Note that the inverse ob jacobian will possibly be
                % updated inside newton, so it mus be a give through
                % parameter.
                [djInv,obj.y(:,k)] = timestep(obj,k-1,djInv,dtime);
            end   
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf(['PDE:BDF2:\nElapsed time ',...
                    num2str(toc),...
                    ' sec using '...
                    num2str(n_jac),' Jacobians.\n']);   
            end
            % Definition of LOCAL Function 
            % This function is NOT the same as in
            % solveEULERI
            function [djInv,y] = timestep(obj,k,djInv,dtime)               
                StepFailed = false;                 
                y = obj.y(:,k); 
                t = obj.time(k+1);
                dy = -obj.df(t,y);
                for k1 = 1:10                  
                    d = -djInv*dy;  
                    [~,id] = lastwarn;
                    if strcmp(id,'MATLAB:nearlySingularMatrix')|| ...
                            strcmp(id,'MATLAB:singularMatrix')
                        StepFailed = true;               
                    end
                    if norm(d,Inf) < 1e-6
                        StepFailed = false;                        
                        break
                    end                
                    y = y + d; 
                    dy = obj.D*((2*dtime(k)+dtime(k-1))/...
                                (dtime(k)^2+dtime(k)*dtime(k-1))*y-...
                                ( dtime(k)+dtime(k-1))/...
                                (dtime(k)*dtime(k-1))*obj.y(:,k)+...
                                ( dtime(k))/...
                                (dtime(k)*dtime(k-1)+dtime(k-1)^2)*obj.y(:,k-1))-...
                                obj.df(obj.time(k),y); 
                    %  Only in case of slow convergence
                    % we compute a new Jacobian,
                    % at 10th 15 and 20th
                    % iteration.
                    switch k1
                        case {10 15 20}
                            djInv = ((2*dtime(k)+dtime(k-1))/...
                                (dtime(k)^2+dtime(k)*dtime(k-1))*obj.D...
                                -obj.jacobian(obj.time(k),obj.y(:,k)))\...
                                speye(size(obj.D));                           
                            n_jac = n_jac+1;
                        otherwise
                            % do nothing                     
                    end                
                end
                if StepFailed                     
                    MException('newton:stepFailed',...                
                        ['Unable to find a solution: ',...
                        lastwarn]).throwAsCaller;                     
                end
                if k1 == MaxNewtonSteps
                    warning('PDE:NewtonSteps',...
                        'Maximal number of Newton Steps reached.')
                end             
            end
        end

        function bdf2sc(obj)
            %bdf2sc (BDF(2) with Stepsize Control)
            % Integrates PDE in time by (implizit) BDF-2 scheme.             
            % A stepsize control is implemented
            %  
            % Usuage:
            % >>obj.bdf2sc()
            % (C) 2013 by Uwe Prüfert
           
            % "local global" focus for jacobian and Step counters
            MaxNewtonSteps  = 50;
            n_jac = 0; 
            n_lu = 0;
            n_failed = 0;
            
            atol = 1e-6;
            rtol = 1e-3;
            
            if isempty(obj.y)||~obj.initialized
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end 
            if isempty(obj.D)
                obj.D = obj.mass;
                warning('PDE:NOTIMEMASS',...
                    'No mass for time derivative given. Use M(d) with d = 1.')
            end
            switch length(obj.time)
                case 0
                     MException(obj.wrongInputFormatID,...
                         [obj.wrongInputFormatStr,...
                         ' obj.time is empty.']).throwAsCaller
                case 1                     
                    obj.time = linspace(0,obj.time,100);
                case 2
                    obj.time = linspace(obj.time(1),obj.time(2),100);
                otherwise
                    % obj.time is assumed as the  vector of timesteps
            end
            
            % First, we will start with a few implicite
            % EULER steps
            % We apply EULERI on the time interval [t(0),t(1)] 
            
            % We must save the original time interval
            
            tmax = max(obj.time);
            
            obj.time = linspace(obj.time(1),obj.time(3),9);
            fprintf('BDF2SC: start with EULER Step\n')
            tic
            obj.euleri;
            
            % extract y0 and the results of euleri
            % at y1 and y2. We need three points
            % to performe the initial predictor
            % step. This is different to the
            % uncontrolled bdf scheme.
            obj.y = obj.y(:,[1,5,9]); 
            obj.time =  obj.time([1,5,9]); 
            % initial step size is dtime
            dtime =  obj.time(2)-obj.time(1);
            % the first steps are chosen equidistant
            % dtime(3) is the guess for the
            % controller!!!
            dtime = dtime(ones(1,3));
            
            % to speed up, we compute here the Jacobian.
            % Inside the loop we use this matrix and we only update it if
            % Newton method converges "slowly".
            % For the newton loop, we set y0 = y(2).
            
            % we store jacobian at y  in J
            
            J = obj.jacobian(obj.time(3),obj.y(:,3));             
            
            % initialized computation of an
            % L-U-decomposition     
            [L,U] = lu((2*dtime(3)+dtime(2))/...
                (dtime(3)^2+dtime(3)*dtime(2))*obj.D-J); 
            
            n_lu = n_lu + 1;
            n_jac = n_jac + 1;
            
            
             
            
            % initialize t by the last time
            % obtained by euler scheme
            t = obj.time(end);
            % we have already three solutions
            k = 3; 
            while t < tmax
                % The time loop, call Newtons methods  to solve the non 
                % linear problem
                % Note that the  jacobian will possibly be
                % updated inside timestep, so it must be a give through
                % parameter.
                
                [L,U,J,obj.y(:,k+1),dtime,t] = timestep(obj,k,L,U,J,dtime,t);               
                obj.time(k+1) = t;
                k = k+1;
            end   
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf(['PDE:solveBDF2:\nElapsed time is ',...
                    num2str(toc),...
                    ' sec.\n',...
                    num2str(n_jac),' Jacobians  and ',num2str(n_lu),...
                    ' LU decompositions were computed for ',num2str(length(obj.time))...
                    ,' time steps.\n',...
                    num2str(n_failed),' steps failed the error bound.\n']);   
            end
            
            
            
            % Definition of LOCAL Function 
            % This function is NOT the same as in
            %  EULERI
            function [L,U,J,y,dtime,t] = timestep(obj,k,L,U,J,dtime,t)  
                JI = J;
                StepFailed = false;             
                while dtime(k) >1e-16
                    y = obj.y(:,k);                     
                    dy = obj.D*((2*dtime(k)+dtime(k-1))/...
                        (dtime(k)^2+dtime(k)*dtime(k-1))*y-...
                        ( dtime(k)+dtime(k-1))/...
                        (dtime(k)*dtime(k-1))*obj.y(:,k)+...
                        ( dtime(k))/...
                        (dtime(k)*dtime(k-1)+dtime(k-1)^2)*obj.y(:,k-1))-...
                        obj.df(1,y); 
                    nd_old = Inf;
                    for k1 = 1:MaxNewtonSteps
                        d = -U\(L\dy);                        
                        [~,id] = lastwarn;
                        if strcmp(id,'MATLAB:nearlySingularMatrix')|| ...
                                strcmp(id,'MATLAB:singularMatrix')
                            StepFailed = true;                              
                        end
                        
                        if norm(d,Inf) < atol
                            StepFailed = false; 
                            break
                        end 
                        
                        
                        y = y + d;
                        % Only in case of slow convergence
                        % we compute a new Jacobian.
                         
                        nd_new = norm(d,Inf);
                        
                        if nd_old/nd_new < 1.4
                             JI = obj.jacobian(obj.time(k),y);
                                [L,U] = lu((2*dtime(k)+dtime(k-1))/...
                                    (dtime(k)^2+dtime(k)*dtime(k-1))*obj.D-JI);
                                n_lu = n_lu + 1;
                                n_jac = n_jac+1; 
                        end 
                        nd_old = nd_new;
                        
                        dy = obj.D*((2*dtime(k)+dtime(k-1))/...
                            (dtime(k)^2+dtime(k)*dtime(k-1))*y-...
                            ( dtime(k)+dtime(k-1))/...
                            (dtime(k)*dtime(k-1))*obj.y(:,k)+...
                            ( dtime(k))/...
                            (dtime(k)*dtime(k-1)+dtime(k-1)^2)*obj.y(:,k-1))-...
                            obj.df(1,y); 
                        
                           
                    end
                   
                    if StepFailed                     
                        MException('newton:stepFailed',...                
                            ['Unable to find a solution: ',...
                            lastwarn]).throwAsCaller;                     
                    end
                    if k1 == MaxNewtonSteps
                        warning('PDE:NewtonSteps',...
                            'Maximal number of Newton Steps reached.');  
                        
                    end
                    
                    % We compute the predictor yP

                    yp = sum(dtime(k-1:k))*sum(dtime(k-2:k))...
                        /dtime(k-1)/sum(dtime(k-2:k-1))*obj.y(:,k)-...
                        dtime(k)*sum(dtime(k-2:k))...
                        /dtime(k-1)/dtime(k-2)*obj.y(:,k-1)+...
                        dtime(k)*sum(dtime(k-1:k))...
                        /dtime(k-2)/sum(dtime(k-2:k-1))*obj.y(:,k-2);  
                    
                    err = norm(y-yp)/...
                        (rtol*norm(y)+atol); 
                    
                    if err<=1
                        % I err  < 1, we can
                        % (i) take y as y(k+1) and update 
                        % t(k+1) = t(k) + dtime(k)
                        % and
                        % (ii) change the stepsize if dt_new differs
                        % from dt_old significantly.                       
                        t = t+dtime(k);                                             
                        dtime(k+1) = dtime(k) * ...
                            min(1+sqrt(2),max(0.3,0.8/(err^3)));
                        
                        if dtime(k+1)/dtime(k) < 1.6 
                            %  We keep dtime(k)
                            dtime(k+1) = dtime(k);                                                  
                        else
                            % otherwise me must update Jacobian but only
                            % the non linear part. Here JI comes into the game         
                            [L,U] = lu((2*dtime(k)+dtime(k-1))/...                                                                                  
                                (dtime(k)^2+dtime(k)*dtime(k-1))*obj.D-JI); 
                            n_lu = n_lu + 1;  
                        end
                        % We can now update J = JI, worst case we have
                        % J = J;-)
                        J = JI;
                        return                       
                    else
                        % Update the Jacobian  because dt has been
                        % changed but only the  linear part...
                        n_failed = n_failed + 1;  
                        dtime(k) = dtime(k) * min(1+sqrt(2),max(0.2,0.9/(err^3)));  
                        % We use here the OLD J...
                        [L,U] = lu((2*dtime(k)+dtime(k-1))/...
                            (dtime(k)^2+dtime(k)*dtime(k-1))*obj.D-J);
                        n_lu = n_lu + 1;                                               
                    end                    
                end
            end
        end           

        function eulerilin(obj)
            %EULERILIN - Integrates LINEAR PDE
            % Integrates LINEAR PDE in time by implizit euler scheme.             
            % An fixed (non adaptive!) time step is used, given by obj.time. 
            % If obj.time contains only t0 te, stepsize is (te-t0)/100 
            % To solve the linear system, pcg is used when Matrix is
            % symmetric otherwise qmr or - if it fails - Gauss-LU will be used.
            % Usuage: Protected method. Use it in your own methods or
            % call it by
            % >>obj.solve('EULERILIN')
            % (C) 2013 by Uwe Prüfert
            tic
             
            maxit = 200;           
            if isempty(obj.y)||~obj.initialized
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end 
            if isempty(obj.D)
                obj.D = obj.mass;
                warning('PDE:solveEULERILIN:NOTIMEMASS',...
                    'No mass for time derivative given. Use M(d) with d = 1.')
            end
            switch length(obj.time)
                case 0
                     MException(obj.wrongInputFormatID,...
                         [obj.wrongInputFormatStr,...
                         ' obj.time is empty.']).throwAsCaller
                case 1                     
                    obj.time = linspace(0,obj.time,100);
                case 2
                    obj.time = linspace(obj.time(1),obj.time(2),100);
                otherwise
                    % obj.time is assumed as the  vector of timesteps
            end           
                                    
            % constant matrix
            % can be a overwritten method...so don't care
            J = sparse((obj.jacobian(obj.time(1),...
                zeros(size(obj.y(:,1),1),1))));            
            % select the inhomogenity (source, boundary conditions 
            % from df...
           
            n = zeros(size(obj.y(:,1),1),1);
            % try symmetric case: pcg
            for k = 2:length(obj.time)
                % time loop,
                dt = obj.time(k)-obj.time(k-1);                     
                % use pcg
                % Algorithm:
                % Note that we here only use the "constant" part from
                % Ay+b by setting y=0 in the call obj.df(t,0).
                % Eventually we want to compute
                % y_k = (D-dt*(K+...+Q))\(D*y_k-1+dt*(R+G+F))
                % Really, it is wrong for non-linear Eqns.
                try
                    L = ichol(sparse(obj.D-dt*J));
                    [obj.y(:,k),flag] = ...
                        pcg(obj.D-dt*J,obj.D*obj.y(:,k-1)+...
                        dt*obj.df(obj.time(k),n),...
                        1e-8,maxit,L,L',obj.y(:,k-1));   
                    
                    
                catch ME
                    % non symmetric case: try qmr
                    % decomposition. 
                    fprintf(['pcg failed. Message was:',ME.message,'\n ']);
                     
                    
                    try
                        [L,U] = ilu(obj.D-dt*J,struct('type','nofill'));
                        [obj.y(:,k),flag] = ...
                            bicgstab(obj.D-dt*J,obj.D*obj.y(:,k-1)+...
                            dt*obj.df(obj.time(k),n),...
                            1e-8,maxit,L,U,obj.y(:,k-1));  
                    catch ME                        
                    % non symmetric case: L-U
                     % decomposition.
                         
                        fprintf(['bicgstab failed. Message was:',ME.message,'\n ']);
                        obj.y(:,k) = (obj.D-dt*J)\(obj.D*obj.y(:,k-1)+...
                            dt*obj.df(obj.time(k),n));
                        flag = 5; 
                    end
                
                end
                switch flag   
                    case 0
                            % everything is fine;-)
                    case 1  
                        fprintf('pcg iterated MAXIT times but did not converge.\n');
                    case 2  
                        fprintf('preconditioner M was ill-conditioned.\n');
                    case 3  
                        fprintf('pcg stagnated (two consecutive iterates were the same.\n');
                    case 4  
                        fprintf(['one of the scalar quantities calculated ',...
                                    'during pcg became too',...
                                    'small or too large to continue computing.\n']);
                    case 5
                        fprintf(['Matrix not symmetric.',...
                            ' Use mldivide operator to solve linear system.\n']);
                    otherwise
                        fprintf('Unknown PCG-flag???\n');
                end
            end 
                 
            % give back the information provided by pcg
            if strcmp(obj.odeOptions.Stats,'on')
                fprintf(['EULERILIN:\nElapsed time ',...               
                    num2str(toc),...
                    ' sec .\n']);   
            end
        end              
        
        function steadystate(obj,te)
            %INITIAL Computes a guess for a initial solution for solveStationary .
            % The result is a vector y(:,te), the time history of the 
            % solution is lost in this way. However, obj.time is now [0, te],
            % where te is the time when obj.eventfun breaks the integration.         
            % solve d dy/dt -f(t,y) = 0 in t = [0,te]
            % te and d are optional.
            % Standard values are d = 1 and te = 1e5. 
            % To stopp the integration, a eventfunction obj.eventfun is
            % used. The standard event to stopp is |df(t,x,y(t,x))|<1e-6.
            % Overwrite it if you wish a different event. The time
            % integration is done by ode15s with standard options.
            % Usuage:
            % >>obj.initial()
            % >>obj.initial(te)
            % (c) 2013 by Uwe Prüfert
            
            switch nargin
                case 1
                    obj.time = [0,1e5];
                     
                case 2
                    obj.time = [0,te];
                otherwise
                    MException(obj.wrongNumberArgumentsID,...
                        [obj.wrongNumberArgumentsStr,...
                        'It must be zero or one.']).thowAsCaller;
            end            
            if isempty(obj.y)||~obj.initialized
                obj.notInitialized.throwAsCaller;
            elseif isvector(obj.y)
                obj.y = obj.y(:);
            else
                obj.y = obj.y(:,end);                            
            end  
            if isempty(obj.D)
                obj.D = obj.mass;
                fprintf(...
                    'No mass for time derivative given. Use Md with d = 1.s\n')
            end
            obj.odeOptions = odeset(...
                obj.odeOptions,...
                    'Mass',obj.D,...
                    'RelTol',1e-6,...
                    'Stats','on',...
                    'Jacobian',@obj.jacobian,...
                    'Events',@obj.eventfun);
            if ~isempty(obj.pattern)
                obj.odeOptions = odeset(obj.odeOptions,...                    
                    'JPattern',obj.pattern);                                
            end
             
            [obj.time,obj.y] = ode15s(@obj.df,[obj.time(1),...
                obj.time(end)/2,obj.time(end)],obj.y,obj.odeOptions);
            obj.y = obj.y(end,:)';  
            obj.time = obj.time(end);  
            fprintf(['System in steady state at time ',...
                num2str(obj.time), '\n']);
        end
        
        function stationary(obj)
            %solveStationary Solves the (non linear) stationary problem. 
            % y can be a sol.y(:,end) from a transient
            % solution process or an arbitrary vector, e.g. y0...
            % dY is a function handles, the same 
            % function as used for ODE15s.
            % The method used a semi Newton method,(without stepsize control) 
            % where the Jacobian is given by obj-jacobian. 
            % The maximum of Newton iterations is 50.
            % If the jacobian is singular or the maximal number of iterations 
            % is reached,the method breaks with an error.
            % If this method fails, especially in advanced cases, feel free
            % to overwrite it. However, to obtain a proper initial guess
            % for y, try a time solver like EULERI
            % Usuage:
            % >>obj.stationary()
            % (c) 2013 by Uwe Prüfert
            if ~obj.initialized
                obj.notInitialized.throwAsCaller
            end
            tic;
            StepFailed = true;  
            d_old = Inf;          
                   
            if isempty(obj.time)
                t = 0;
            else
                t = obj.time(end);
            end
             
            if isempty(obj.y)
                % Need it for compute the jacobian
                dimProb = max(size(obj.A,1),size(obj.K,1));
                obj.y = zeros(dimProb,1);
            else
                obj.y = obj.y(:,end);
            end 
            dy = obj.df(t,obj.y);
            dj = -obj.jacobian(t,obj.y);           
           
            n_jac = 1;
            for k1 = 1:50   
                
                d = dj\dy;  
                [~,id] = lastwarn;
                if strcmp(id,'MATLAB:nearlySingularMatrix')|| ...
                        strcmp(id,'MATLAB:singularMatrix')                    
                    StepFailed = true;               
                end
                if norm(d,Inf) < obj.solverOptions.solverTol  
                    fprintf(['Newton loop breaks at ',...
                        num2str(k1),'-th iteration.\n']); 
                    StepFailed = false;
                    break
                end
                obj.y = obj.y + d;
                dy = obj.df(t,obj.y); 
                d_monitor = norm(d-d_old);
                d_old = d;
                % for large different d's compute new Jacobian
                if d_monitor > .1 
                    dj = -obj.jacobian(t,obj.y);
                    n_jac = n_jac+1;
                end                
            end
            if StepFailed
                
                ME = MException('stationary:stepFailed',...                
                    'Unable to find a stationary solution. ');                
                ME.throwAsCaller;
            else
                fprintf('Success: I found a stationary solution ')
                fprintf(['within ',num2str(toc),' sec.\n']);
                fprintf(['Number of Jacobians to compute ',...
                    num2str(n_jac),'.\n\n']);            
            end
        end
        
        function linear(obj) 
            %linear solves dF = 0 by computing y = JdF\dF Hence it
            % solves a stationary, linear PDE  problem.
            % In contrast to all other solvers, obj.y can be empty. It will
            % be overwritten by zeros and - in the case of successful
            % solving - by the solution of the linear problem.
            % The user should provide obj.jacobian - in this case a
            % function gives back a constant matrix - to speed up solution
            % proccess. Otherwise, obj.jacobian will compute this matrix in
            % a more time consuming way.
            % Usuage: It is a protected method so you can use it in
            % your own methods or call it by using
            % >>obj.solve('LINEAR')
            % (c) 2013 by Uwe Prüfert
            maxit = 200;
            if ~obj.initialized
                obj.notInitialized.throw
            end
            tic
            if isempty(obj.y)
                % Need it for compute the jacobian
                dimProb = max(size(obj.A,1),size(obj.K,1));
                obj.y = zeros(dimProb,1);
            else
                obj.y = obj.y(:);
            end 
            if isempty(obj.time)
                t = 0;
            else
                t = obj.time(end);
            end
            try               
                L = ichol(-sparse(obj.jacobian(t,obj.y)));            
                [obj.y,flag] = ...
                    pcg(sparse(-obj.jacobian(t,obj.y)),obj.df(t,obj.y),...
                    1e-8,maxit,L,L');                
            catch ME
                % non symmetric case: L-U
                 % decomposition.
                fprintf(['iterative solver failed: ',ME.message,...
                    '\nTry Gauss LU\n']);
                obj.y = -obj.jacobian(t,obj.y)\obj.df(t,obj.y);
                flag = 5; 
            end
            switch flag   
                case 0
                     fprintf('  pcg converged.\n');   % everything is fine;-)
                case 1  
                    fprintf('pcg iterated MAXIT times but did not converge.\n');
                case 2  
                    fprintf('preconditioner M was ill-conditioned.\n');
                case 3  
                    fprintf('pcg stagnated (two consecutive iterates were the same.\n');
                case 4  
                    fprintf(['one of the scalar quantities calculated ',...
                                'during pcg became too',...
                                'small or too large to continue computing.\n']);
                case 5
                    fprintf(['Matrix not symmetric.',...
                        ' Use mldivide operator to solve linear system.\n']);
                     [~,id] = lastwarn;
                     if strcmp(id,'MATLAB:nearlySingularMatrix')|| ...
                                     strcmp(id,'MATLAB:singularMatrix')
                        fprintf(['Solved linear problem with warnings ',...
                                id,'\n']);
                     end
                otherwise
                    fprintf('Unknown PCG-flag???\n');
            end
            fprintf(['Spend ',...
                    num2str(toc),' seconds for solving pde\n'])
             
        end
        
        function linearGauss(obj)
            %linearGauss solves dF = 0 by computing y = JdF\dF Hence it
            % solves a stationary, linear PDE  problem similar to LINEAR 
            % but with Gauss LU decomposition.
            % In contrast to all other solvers, obj.y can be empty. It will
            % be overwritten by zeros and - in the case of successful
            % solving - by the solution of the linear problem.
            % The user should provide obj.jacobian - in this case a
            % function gives back a constant matrix - to speed up solution
            % proccess. Otherwise, obj.jacobian will compute this matrix in
            % a more time consuming way.
            % Usuage: It is a protected method so you can use it in
            % your own methods or call it by using
            % >>obj.solve('LINEAR')
            % (c) 2014 by Uwe Prüfert
         
            if ~obj.initialized
                obj.notInitialized.throw
            end
            tic
            if isempty(obj.y)
                % Need it for compute the jacobian
                dimProb = max(size(obj.A,1),size(obj.K,1));
                obj.y = zeros(dimProb,1);
            else
                obj.y = obj.y(:);
            end 
            if isempty(obj.time)
                t = 0;
            else
                t = obj.time(end);
            end
            try               
                obj.y = -obj.jacobian(t,obj.y)\obj.df(t,obj.y);
            catch ME
                ME.thowAsCaller                
            end
            fprintf(['Spend ',...
               num2str(toc),' seconds for solving pde\n'])
             
        end
        
    end
    
    methods(Access=protected)
        %df RHS of equation d*dy/dt = df(t,y)
        % eg. df(t,y) = -(K+M+Q)*y+(G+F) for a heat transfer problem with
        % Robin BCs.
        % The method df is here declarated as ABSTRACT.
        % Hence, the df method must be implemented in every child class!
        % THIS IS A DEMAND ON USER INPUT FOR THE CHILD CLASSES
        % Signature: 
        %    dY = df(obj,t,y)  
        % In child classes, it must have Access = protected.
        % (c) 2013 by Uwe Prüfert
        dY = df(obj,t,y)  
    end
    
    methods(Static)
    end
    
    methods(Hidden,Access=public)        
        function addlistener(obj)
            obj.notAllowed.throwAsCaller;
        end
        function eq(obj,~)
            obj.notAllowed.throwAsCaller;
        end 
        function ne(obj,~)
            obj.notAllowed.throwAsCaller;
        end
        function lt(obj,~)
            obj.notAllowed.throwAsCaller;
        end 
        function le(obj,~)
            obj.notAllowed.throwAsCaller;
        end 
        function ge(obj,~)
            obj.notAllowed.throwAsCaller;
        end 
        function gt(obj,~)
            obj.notAllowed.throwAsCaller;
        end 
        function findobj(obj,~)
            obj.notAllowed.throwAsCaller;
        end
        function findprop(obj,~) 
            obj.notAllowed.throwAsCaller;
        end
        function notify(obj,~)
            obj.notAllowed.throwAsCaller;
        end
    end
end