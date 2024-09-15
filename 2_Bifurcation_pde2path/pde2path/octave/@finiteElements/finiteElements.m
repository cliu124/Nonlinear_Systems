classdef  finiteElements < handle
    % finiteElements
    % Abstract class to build Finite Elements    
    % It defines only universal error Objects.
    
    properties(Constant)
        wrongNumberInputs = MException('FINITEELEMENTS:WRONGNUMBERINPUTS',...
                        'The number of input arguments is wrong.');
        wrongNumberOutputs = MException('FINITEELEMENTS:WRONGNUMBERINPUTS',...
                        'The number of output arguments is wrong.');         
        wrongFormat = MException('FINITEELEMENTS:WRONGFORMAT',...
                                    'reading geometry data failed, maybe incorrect format.');
        wrongClass = MException('FINITEELEMENTS:WRONGCLASS',...
                                    'Wrong argument class.');
        wrongSize = MException('FINITEELEMENTS:WRONGSIZE',...
                                    'Wrong sized arguments.');                        
    end  
  
    %properties(Constant); idx; end
    
    methods(Access=public)
        function M = sourceTermMatrix(obj,grid) 
    %%
    % * sourceTermMatrix IN:self,gridd OUT: double
    %
    % Computes the matrix M such that F = M*f,
    % where f ist a vector of length #Elements which contrains the
    % values of the source f in the centers of each element. This
    % is usefull if you want to calculate the source from a given
    % vector, e.g.  when solving coupled systems where the solution
    % of one equation is the source for another equation.
    %
    % Call:
    %
    %        M = finiteElements.sourceTermMatrix(gridd)

    % (c) 2016 Uwe Pruefert.

    idx1 = reshape(grid.t(1:grid.nPointsInElements,:),...
        1,grid.nElements*grid.nPointsInElements);
    idx2 = reshape(repmat(1:grid.nElements,grid.nPointsInElements,1),...
        1,grid.nElements*grid.nPointsInElements); 

    M = sparse(idx2,idx1,reshape(obj.F*obj.makeJ(grid),1,grid.nPointsInElements*grid.nElements),grid.nElements,grid.nPoints)';
        end
        
    function errorPerElement = errorInd(obj,gridObj,y,c,a,f,alpha,beta,m)
        %size(c), size(a), size(f)
           switch nargin
                case 6
                    alpha = 1;   beta = 1;   m = 2;
                case 8
                    m = 2; 
                case 9
                    % okay
                otherwise
                    gridObj.wrongNumberInputs.throwAsCaller;
           end
            [sideLength,~]= gridObj.sideLengthAndArea();
            normFminusau = obj.localErrorL2(gridObj,y,a,f); 
            
            % multiply by triangle's longes side 
            normFminusau = normFminusau.*max(sideLength).^m; 
            % fluxes through edges, static and abstract method
            %try
            fluxThroughElementEdges = obj.fluxThroughEdges(gridObj,y,c);           
           
            
            % computeFluxJumps static and abstract method
            jumps = obj.fluxJumps(gridObj,fluxThroughElementEdges,m);
            %catch jumps=0; end % HU
            % okay, all methods are abstract and static and should work for 1D--3D
            % -------------------------------------------------------------           
          
            errorPerElement = alpha*normFminusau + beta*sqrt(0.5*jumps);
    end
    end
    methods(Static,Access = public)
        % Functions for discretization error handling: adaptive error
        % indicator etc. NOTE: Tricky! These methods are static. The argument
        % gridObbj may differ from element-type to element-type.
        % fluxThrougElementEdges= fluxThroughEdges(gridObj,u,c) 
        % jumps = fluxJumps(gridObj,fluxThroughElementEdges,order)
    end 
     methods(Static,Access = public) 
        % localErrorL2 is made for local error measurment. It may be
        % overwritten in higher order FE classes. However, for linear FEs
        % it works for every dimension
        function localL2 = localErrorL2(obj,u,a,f)
            %localErrorL2  L2 norm of f-a*u
            % Evaluates f-a*u in the center of Element and multiplies with
            % area. Result is per element. This should work for all linear
            % P1 elements in 1D--3D. May be overwritten in P2 etc. classes.
            % error = obj.localErrorL2(u,a,f).
            %
            %Note that u, a and f  must be vectors or scalars, not function
            %handles. 
            % (c) 2015 Uwe Prüfert

            % Changelog May 4th 2105: Make it dimension independent

            % this must be in the gridxD classes. In 1D and 3D may area be
            % misleading: It is the volume of the lement and length is the
            % length or area of the edge or boundary 
            [~,area] = obj.sideLengthAndArea;
            %size(f), size(u), size(a), pause % HU 
            if isscalar(f); f=ones(1,obj.nElements)*f; 
            elseif min(size(f))==1&&max(size(f))==obj.nPoints
                f = obj.point2Center(f); f = f(:)'; 
            elseif min(size(f))==1&&max(size(f))==obj.nElements; f = f(:)';                
            else MException(obj.wrongInputFormatID,...
                    [obj.wrongInputFormatStr,...
                    ' must be vector of lenght np or ne,  or must be a  scalar.']).throwAsCaller;
            end
            if isscalar(a); a=ones(1,obj.nElements)*a; % HU, was f
            elseif min(size(a))==1&&max(size(a))==obj.nPoints
                a = obj.point2Center(a); a = a(:)';   
            elseif min(size(a))==1&&max(size(a))==obj.nElements
                a = a(:)'; 
            else MException(obj.wrongInputFormatID,...
                    [obj.wrongInputFormatStr,...
                    ' must be vector of lenght np or ne,  or must be a  scalar.']).throwAsCaller;
            end 
    if min(size(u))==1 && max(size(u))==obj.nPoints; u=obj.point2Center(u); u=u(:)'; end % HU 
            if size(u,1)>size(u,2); u=u'; end 
           % 13, size(f), size(u), size(a), pause 
            localL2 = sqrt((f-a.*u).^2.*area); 
            
        end
    end
    methods(Access = public)   
        % Assembling of all Matrices valid for all ...        
        function [varargout] = assema(obj,gridObj,cf,af,ff)
        % [K,M,F] = fem.assema(grid,c,a,f)
        % [S,F] = fem.assema(grid,c,a,f)
        %  The assemble method for finite elements.
        %  It assembles the Stiffness matrix K, the Mass Matrix M and the RHS 
        %  vector F in sparse format. There is a two argument call where S = K+M.  
        % (c) 2013 by Uwe Prüfert
            % specialized methods dep. on element type
            
           
            [K,M,F] = obj.createMatrixEntries(gridObj,cf,af,ff); 
            [idx0,idx1,idx2] =  obj.makeIndex(gridObj.t(obj.idx,:),gridObj.nElements); 
            % do the sparse magic...  
            switch nargout
                case 2  
                    varargout{1} = sparse(idx1,idx2,K+M,gridObj.nPoints,gridObj.nPoints); 
                    varargout{2} = sparse(idx0,1,F,gridObj.nPoints,1);
                case 3  
                    varargout{1} = sparse(idx1,idx2,K,gridObj.nPoints,gridObj.nPoints);  
                    varargout{2} = sparse(idx1,idx2,M,gridObj.nPoints,gridObj.nPoints);
                    varargout{3} = sparse(idx0,1,F,gridObj.nPoints,1);
                otherwise
                    throw(obj.wrongNumberOutputs)
            end         
        end
        
        function B = convection(obj,gridObj,b)
        % convection method of class finiteElements.
        % Universal convection assembling
        % B = fe.convection(grid,bvec) 
        % (C) 2013 by Uwe Prüfert.            
%             idx = gridObj.t(obj.idx,:);   
            % makeIndex method call, this method is in this class abstract! 
            
            [~,idx1,idx2] = obj.makeIndex(gridObj.t(obj.idx,:),gridObj.nElements);
            val = obj.createConvectionEntries(gridObj,b);
            B = sparse(idx1,idx2,val,gridObj.nPoints,gridObj.nPoints);    
        end
        
        
        function B = sparsityPattern(obj,gridObj,varargin)
            % sparsityPattern computes the sparsity pattern of the system
            % matrices
            % Usuage:
            % B = obj.sparsityPattern(gridObject)
            % B = obj.sparsityPattern(gridObject,A)
            % B = obj.sparsityPattern(gridObject,K,M,C,Q,H)
            % A can be every matrix, 
            % K,M,C,Q,H  should be the matrices from a call of obj.assema,
            % obj.assemb and obj.convection
            % The empty call computes all relevant matrices internally 
            
            switch length(varargin)
                case 0  % no matrices given
                    [K,M,~] = obj.assema(gridObj,'1','1','1');
                    [Q,~,H,~] = obj.assemb(gridObj);
                    C = obj.convection(gridObj,ones(size(gridObj.p(:,1))));
                    [indx1,indx2,s] = find(K+M+C+(H'*H)+Q);
                    B = sparse(indx1,indx2,s~=0);                    
                case 1  % whole matrix
                    [indx1,indx2,s] = find(varargin{1});
                    B = sparse(indx1,indx2,s~=0);                   
                case 5  % all relevant matrices in order K,M,C,Q,H
                    [indx1,indx2,s] = find(varargin{1}+...
                                           varargin{2}+...
                                           varargin{3}+...
                                           varargin{4}+...
                                           varargin{5}'*varargin{5});
                    B = sparse(indx1,indx2,s~=0); 
                otherwise
                    obj.wrongNumberInputs.throwAsCaller;
            end
        end
    end
    
    methods(Static,Access=public)
        function L = stiffSpring(M)
            % computes a guess for the stiff-spring coefficient
            % L = obj.stiffSpring(M)
            % M should be the Systems Matrix, e.g. M = K+M+C
            % (C) 2013 by Uwe Prüfert
            L = 1e3*norm(M,1);
        end
    end
    
    % declarations
    methods(Access = public)   
        % Assembling of  Matrices abstract definition
        % to be specialized in subclasses
                  
     %   [K,M,F] = createMatrixEntries(obj,gridObj,cf,af,ff);
     %   val = createConvectionEntries(obj,gridObj,b);
     %   [idx0,idx1,idx2] = makeIndex(idx,nElements);
    end 
    
    methods(Static,Access = public)   
        % Assembling of  Matrices abstract definition
        % to be specialized in subclasses
       % [Q,G,H,R] = assemb(gridObj);                  
        
    end 
end

