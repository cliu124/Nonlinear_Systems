% finiteElements2D
classdef finiteElements2D < finiteElements
    properties(Access = public)
        % boundaryElements finiteElements
        % Defines the boundary elments for assembling the boundary
        % condtion. 
    end
    
    properties(Access = public)
        % The order of point indeces in the e field. Trivialy [1,2] for P1 
        % elements but importand for extended meshes.
       % boundaryIndex double;
    end
        
    methods(Access = public)     
        function [K,M,F] = createMatrixEntries(obj,gridObj,cf,af,ff)
            % the sizes of elementary matrices, varies from element
            % type to element type...
            % The space needed by sparse methos. Note that S1 etc are
            % stored as vectors. It length is 9 in the case of linear
            % Elements on triangles but it can vary 
            sizeVector = sqrt(size(obj.S1,1));
            sizeMatrix = size(obj.S1,1);
            % compute the coefficients. The result is always a k x nElelemts
            % matrix where k = 3 for c, k = 1 for a and f.
            [cval,aval,fval] = obj.aCoefficients(gridObj,cf,af,ff);
            switch size(cval,1) 
                case 2
                    % We obtain the Matrices from the exact integration multiplied
                    % by the grid-dependent but for fix grid constant
                    % transformations.
                    [K1,~,K2] = obj.getConstantPartOfStiffnessMatrix(gridObj);           
                    % Now we have to multiply the Ks by c(1) and c(2).
                    % For the mixed derivatives, this code must be EXTENDED.
                    K = reshape(...
                            K1*sparse(1:gridObj.nElements,1:gridObj.nElements,...
                                cval(1,:))+ K2*sparse(1:gridObj.nElements,...
                                1:gridObj.nElements,...
                                cval(2,:)),1,gridObj.nElements*sizeMatrix);

                case 4                    
                    [K1,K12,K2] = obj.getConstantPartOfStiffnessMatrix(gridObj);           
                     K = reshape(...
                        K1*sparse(1:gridObj.nElements,1:gridObj.nElements,...
                        cval(1,:))+ K2*sparse(1:gridObj.nElements,1:gridObj.nElements,...
                        cval(4,:))+ K12*sparse(1:gridObj.nElements,1:gridObj.nElements,...
                        cval(3,:)+cval(2,:)),1, gridObj.nElements*sizeMatrix);
                otherwise
                    error(['The object cval has ' num2str(size(cval,1)) ' Rows.'...
                        'If Diffusion Coefficient is scalar, cval must have three equal rows. '...
                        'If it is a diagonal matrix, cval must have these 2 diagonalelements as rows. '...
                        'If it is a full matrix, cval must have 4 rows, which are taken from the matrix columnwise. '...
                        'So if the rows are labeld with c1,c2,c3,c4, then the matrix entries will be'...
                        ' [c1 c3; c2 c4].'...
                        ' Other syntaxes are not supported. For further help, see aCoefficients() in your used grid class.']);
            end
                    
            M = reshape(obj.M*(aval.*obj.makeJ(gridObj)),1,gridObj.nElements*sizeMatrix);             
            F = reshape(obj.F*(obj.makeJ(gridObj).*fval),1,gridObj.nElements*sizeVector);
        end
        
        function val = createConvectionEntries(obj,gridObj,b)
            % the point of the triagles
            p1 = gridObj.p(:,(gridObj.t(1,:)));
            p2 = gridObj.p(:,(gridObj.t(2,:)));
            p3 = gridObj.p(:,(gridObj.t(3,:)));
            % get element corner coordinates
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:);
            bval = obj.convCoefficients(gridObj,b);
            J1 = bval(:,1)'.*y31-bval(:,2)'.*x31; 
            J2 = bval(:,2)'.*x21-bval(:,1)'.*y21;      
            % This works for P1 and P2 etc. elements, really
            val = reshape(obj.C1*J1+obj.C2*J2,1,gridObj.nElements*size(obj.C1,1));             
        end          
    end
    
    methods(Access = public,Hidden)
        function [K1,K12,K2] = getConstantPartOfStiffnessMatrix(obj,gridObj)
            % Needed by createMatrixEntries
            % helper method to compute the former abstract J 
            %                        
            p1 = gridObj.p(:,(gridObj.t(1,:)));
            p2 = gridObj.p(:,(gridObj.t(2,:)));
            p3 = gridObj.p(:,(gridObj.t(3,:))); 
            
            % get element corner coordinates
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:); 
            
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:); 
            
            % now we compute J for every element in one operation              
            J = obj.makeJ(gridObj) ;    
            % Inversion by using the famous 2x2 inversion formula, but with
            % out division by J: A few lines later we compute e.g.
            % xi_x^2*j so we can "correct" this "mistake"  here.
            xi_x  =  y31; 
            eta_x = -y21;
            xi_y = -x31;
            eta_y = x21;  
            
            K1 = obj.S1*((xi_x.^2)./J) + obj.S3*((eta_x.^2)./J)...
                +  obj.S2*((xi_x.*eta_x)./J);
            
            K2 = obj.S1*((xi_y.^2)./J) + obj.S3*((eta_y.^2)./J)...
                +  obj.S2*((xi_y.*eta_y)./J);
            % The following will work only for symmetric matrices.
            K12 = obj.S1*((xi_y.*xi_x)./J) + obj.S3*((eta_y.*eta_x)./J)...
                +  obj.S2*((xi_y.*eta_x+ xi_x.*eta_y )./J/2);   
        end
    end
    
    methods(Static,Access=public)
        function [cval,aval,fval] =  aCoefficients(gridObj,cc,aa,ff)
            [cval,aval,fval] = gridObj.aCoefficientsMpt(cc,aa,ff);
        end       
        
        function bval = convCoefficients(gridObj,b)
            bval = gridObj.convCoefficientsMpt(b);
        end
     
        function J = makeJ(gridObj)   
            % compute Jacobi determinat
            if ~(isa(gridObj,'grid2D')||isa(gridObj,'grid2DR'))
                gridObj.wrongClass.throwAsCaller;
            end
            %heim_neu
            if gridObj.isExtended
                p = gridObj.p(:,1:gridObj.NrPO);
            else
                p = gridObj.p;
            end
            % end heim_neu
            t = gridObj.t;
            p1 = p(:,(t(1,:)));
            p2 = p(:,(t(2,:)));
            p3 = p(:,(t(3,:)));
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:);            
            J = x21.*y31-x31.*y21;               
        end
        
        
        %Implement abstract inheritance from finiteElements class
        function ddncu = fluxThroughEdges(obj,u,c)
            % fluxThroughEdges Fluxes of -div(c grad(u)) 
            % through edges of triangles.
            % ddncu = obj.fluxThroughEdges(grid,u,c)

            
            % Arguments check:
            switch length(c)
                case {1, obj.nElements}
                    % c fits!
                case obj.nPoints
                    c = obj.point2Center(c);
                otherwise
                    % c is miss-dimensioned
                    obj.wrongFormat.throwAsCaller;  
            end
            if length(u)~=obj.nPoints
                obj.wrongFormat.throwAsCaller; 
            end
                
            
            [sideLength, elementArea] = obj.sideLengthAndArea;
            

            dx = zeros(3,obj.nElements); dy = zeros(3,obj.nElements);
            for k = 1:3
                k1 = rem(k ,3)+1;
                k2 = rem(k1,3)+1;
                dx(k,:) = obj.p(1,obj.t(k1,:)) - obj.p(1,obj.t(k2,:));
                dy(k,:) = obj.p(2,obj.t(k1,:)) - obj.p(2,obj.t(k2,:));
            end
            % gradients of basis functions
            g1x=0.5*dy(1,:)./elementArea;
            g2x=0.5*dy(2,:)./elementArea;
            g3x=0.5*dy(3,:)./elementArea;
            g1y=-0.5*dx(1,:)./elementArea;
            g2y=-0.5*dx(2,:)./elementArea;
            g3y=-0.5*dx(3,:)./elementArea;

            % preparations  

            % select points from triangle
            it1=obj.t(1,:);
            it2=obj.t(2,:);
            it3=obj.t(3,:);

            % Compute gradient numerically
            % grad_x
            gradx = u(it1,:).'.*(g1x)+u(it2,:).'.*(g2x)+ ...
                  u(it3,:).'.*(g3x);

            % grad_y  
            grady = u(it1,:).'.*(g1y)+u(it2,:).'.*(g2y)+ ...
                  u(it3,:).'.*(g3y);

            % Compute c * grad u
            nrc=size(c,1);
            cgradx=zeros(1,obj.nElements);
            cgrady=zeros(1,obj.nElements);
            if nrc==1 
                cgradx(1,:)=c.*gradx(1,:);
                cgrady(1,:)=c.*grady(1,:);                
            elseif nrc==2 
                cgradx(1,:)=c(1,:).*gradx(1,:);
                cgrady(1,:)=c(2,:).*grady(1,:); 
            elseif nrc==3 
                cgradx(1,:)=c(1,:).*gradx(1,:)+c(2,:).*grady(1,:);
                cgrady(1,:)=c(2,:).*gradx(1,:)+c(3,:).*grady(1,:);               
            elseif nrc==4 
                cgradx(1,:)=c(1,:).*gradx(1,:)+c(3,:).*grady(1,:);
                cgrady(1,:)=c(2,:).*gradx(1,:)+c(4,:).*grady(1,:);
            else
              MException(obj.wrongInputFormatID,...
                    [obj.wrongInputFormatStr,'Wrong format of c.']).throwAsCaller
            end

            % nhat'*c grad u
            % edge unit normals : outwards positive if the nodes are in
            % anti-clockwise order
            % nhatx =   dy./s;
            % nhaty = - dx./s;
            ddncu = zeros(3,obj.nElements);             
            for k = 1:3
                ddncu(k,:) = (dy(k,:).*cgradx(1,:) - dx(k,:).*cgrady(1,:))./sideLength(k,:);
            end
        end

        function jumps = fluxJumps(obj,fluxThroughElementEdges,m)
            jumps = zeros(1,obj.nElements);
            [sideLength, ~] = obj.sideLengthAndArea;

            % --- this is dimension depending -- here 2D triangle
            interiorEdges=sparse(obj.t([2 3 1],:),obj.t([3 1 2],:),1,obj.nPoints,obj.nPoints);
            % intj+intj.' is 2 if interior edge and 1 if exterior edge

            % --- Interior only
            interiorEdges=round((interiorEdges+interiorEdges.')/3);



            jmps = sparse(obj.t([2 3 1],:),obj.t([3 1 2],:),...
                fluxThroughElementEdges([1 2 3],:),obj.nPoints,obj.nPoints);
            

            jmps = interiorEdges.*abs(jmps + jmps.');
            for l = 1:obj.nElements
                jumps(l) = (sideLength(3,l)^m*abs(jmps(obj.t(1,l),obj.t(2,l))))^2+...
                       (sideLength(1,l)^m*abs(jmps(obj.t(2,l),obj.t(3,l))))^2+...
                       (sideLength(2,l)^m*abs(jmps(obj.t(3,l),obj.t(1,l))))^2;
            end
        end    
    end
    
    methods(Access=public,...
            Static = true)
        

        function [H,R] = assembDataDirichlet(gridObj,arg1,arg2)
            
            switch nargin
                case 1 % No data given: Assume homogenious Dirichlet Problem 
                    arg1 = 1;
                    arg2 = 0;
                case 2 % only R given: Implement h = r
                    arg2 = arg1;
                    arg1 = 1;
                case 3
                    % Nothig to do.
            end
            
            if ~(isa(gridObj,'grid2DR')||isa(gridObj,'grid2D'))
                MException('BILINEAR3D:WRONGARGS',...
                        'The first argument must be a grid2DR object.').throwAsCaller;
            end
            
            % The case of intersecting Dirichlet boundary segments.
            if gridObj.nDirichletBoundaryPoints~=...
                length(gridObj.indexOfDirichletBoundaryPoints)
                fprintf(['\nWarning: Two or more Dirichlet boundary segments have common points.\n',...
                    'In this case, the boundary value function may be not continuous.\n',...
                    'The value that appears first in the arguments \n',...
                    'will be taken, the second will be ignored.\n\n']);
                
                [~,s,~] = unique(gridObj.indexOfDirichletBoundaryPoints,'stable');
                              
                if isscalar(arg1)
                    h = arg1(ones(1,gridObj.nDirichletBoundaryPoints));
                else
                    h = arg1(s);
                end
                
                if isscalar(arg2)
                    R = arg2(ones(gridObj.nDirichletBoundaryPoints,1));
                else
                    R = arg2(s);
                end
                
                H = sparse(1:gridObj.nDirichletBoundaryPoints,...
                        gridObj.indexOfDirichletBoundaryPoints(s),...
                        h,...
                        gridObj.nDirichletBoundaryPoints,...
                        gridObj.nPoints);  
                
                
                R = R(:);
            else
                if isscalar(arg1)
                    h = arg1(ones(1,gridObj.nDirichletBoundaryPoints));
                else
                    h = arg1;
                end
                
                if isscalar(arg2)
                    R = arg2(ones(gridObj.nDirichletBoundaryPoints,1));
                else
                    R = arg2(:);
                end
                H = sparse(1:gridObj.nDirichletBoundaryPoints,...
                        gridObj.indexOfDirichletBoundaryPoints,...
                        h,...
                        gridObj.nDirichletBoundaryPoints,...
                        gridObj.nPoints);  
            end
        end
    end
    
    methods(Access=public)
        function [Q,G] = assembDataRobin(obj,gridObj,arg1,arg2)
            %% assembDataRobin
            %  IN:self, grid2DR,double[,double] OUT:double,double
            %
            % Assembles matrix Q and G from data vectors.
            % Arguments are the grid, the data g (Neumann BC) or q and g.
            % 
            % The data vectors must contain the value of q(x) and g(x)
            % respectivelly in the barycenters of the boundary elements,
            % triangle or sqares, respectivelly. 
            %
            % The order of the data must be the order of elements in obj.e.
            % of the Robin and/or Neumann baundary conditions.
            
            % prepare and test data.
            
            % Since the boundary elements in the boundary  segments are not
            % increasing we resort the data.
            switch nargin
                case 3 % Only g given: Neumann boudarys 
                    g = arg1;    
                    if isscalar(g)
                        g = g(ones(1,gridObj.nRobinBoundaryElements));
                    elseif isvector(g)&&length(g)==gridObj.nRobinBoundaryElements  
                            %  okay
                    else
                        gridObj.wrongFormat.throwAsCaller;
                    end
                    q = zeros(1,gridObj.nRobinBoundaryElements);
                case 4
                    g = arg2;
                    if isscalar(g)
                        g = g(ones(1,gridObj.nRobinBoundaryElements));
                    elseif isvector(g)&&length(g)==gridObj.nRobinBoundaryElements  
                        %  okay 
                    else
                        gridObj.wrongFormat.throwAsCaller;
                    end                     
                    q = arg1;
                    if isscalar(q)
                        q = q(ones(1,gridObj.nRobinBoundaryElements));
                    elseif isvector(q)&&length(q)==gridObj.nRobinBoundaryElements  
                         % okay
                    else
                        gridObj.wrongFormat.throwAsCaller;
                    end
                                      
                otherwise %
                    MException('BILINEAR3D:WRONGNUMBERARGS',...
                        'Wrong number of arguments.').throwAsCaller;
            end
            if ~(isa(gridObj,'grid2DR')||isa(gridObj,'grid2D'))
                MException('finiteElements2D:WRONGARGS',...
                        'The first argument must be a grid2DR object.').throwAsCaller;
            end
                        
            %  index of Faces 
%             indexOfRobinBEs = gridObj.indexOfRobinBoundaryElements;
            
            % x-y-coordinates of the three face points... 
            % p = gridObj.p(1:2,gridObj.e(1:3,k));
            p1 = gridObj.p(1:2, gridObj.e(1,gridObj.indexOfRobinBoundaryElements));
            p2 = gridObj.p(1:2, gridObj.e(2,gridObj.indexOfRobinBoundaryElements));            

            % get element corner coordinates
            x21 = p2(1,:)-p1(1,:);            
            y21 = p2(2,:)-p1(2,:);                 
            
            % Compute Jacobi determinat wrt. boundary points. 
            J = sqrt(x21.^2+y21.^2);          
            
            % Use elementary mass matrix from objboundaryElement class to compute
            % the 1 dimensional integrals. 
            
            dm = length(obj.boundaryElements.M);
            df = length(obj.boundaryElements.F);
            
            Qe = reshape(obj.boundaryElements.M*(J.*q),1,dm*gridObj.nRobinBoundaryElements);
            ge = reshape(obj.boundaryElements.F*(J.*g),1,df*gridObj.nRobinBoundaryElements); 
            
            % rearanging like in assema...
            % Note that we must take case if the mesh is extended.
            if gridObj.isExtended
                 indexOfPoints = gridObj.e([1 2 6],gridObj.indexOfRobinBoundaryElements);
            else
                indexOfPoints = gridObj.e(1:2,gridObj.indexOfRobinBoundaryElements);
            end
 
            % We use again obj.booundaryElements vlass to compute
            % indeces...
            [indx0,indx1,indx2] = obj.boundaryElements.makeIndex(indexOfPoints,gridObj.nRobinBoundaryElements) ;       
            
            % Similar to assema: sparse-magic
            Q = sparse(indx1,indx2,Qe,gridObj.nPoints,gridObj.nPoints); 
            G = sparse(indx0,1,ge,gridObj.nPoints,1);              
        end
        
        function [Q,G,H,R] = assemb(obj,gridObj)            
            % HU: old version
            Q = spalloc( gridObj.nPoints, gridObj.nPoints,gridObj.nPoints);
            H = spalloc( gridObj.nPoints, gridObj.nPoints,gridObj.nPoints);
            G = spalloc( gridObj.nPoints,1, gridObj.nPoints); 
            R = spalloc( gridObj.nPoints,1,  gridObj.nPoints);
                        
            p = gridObj.p; e = gridObj.e; b = gridObj.b;
                           
            for k = 1:gridObj.nEdges
                % length of the egde segment
                l = sqrt((p(1,e(2,k))-p(1,e(1,k)))^2 + ...
                    (p(2,e(2,k))-p(2,e(1,k)))^2);
                
                % simple L1 case, same as in assemb
                % first testing if coefficients are parametrized by
                % arclength
                if is_arclength_parametrized(b(:,e(5,k)))
                    midpoint_of_the_edge = 0.5*( e(3,k)+e(4,k));
                else
                    midpoint_of_the_edge = 0.5*[p(1,e(2,k))+p(1,e(1,k))...
                        ,p(2,e(2,k))+p(2,e(1,k))];
                end
                % where we put the entries...
                indx = e(1:2,k);
                
               [qval,gval,hval,rval] = gridObj.boundCoefficients(midpoint_of_the_edge,b(:,e(5,k)));
             %      [qval,gval,hval,rval] = gridObj.boundCoefficients(e(5,k),b(:,e(5,k)));  % HU 
                
                % Here we postpone the "clever sparse action" and accept a
                % slower code...
                Q(indx, indx) = Q(indx, indx)+...
                    qval*l*reshape(lagrange11D.M,2,2); %#ok<*SPRIX>
                G(indx) = G(indx)+gval*l*lagrange11D.F; 
                H(indx, indx) = H(indx,indx)+hval*diag([1 1],0);
                R(indx) = R(indx)+rval*[1;1];     
            end
              function xy = is_arclength_parametrized(b)             
                % checks if the function in b are arc-length-parametrized  
                % or not by trial-and-error.
                % Idee: If (x,y) is given and b is a function of b(x,y) 
                % or b(x) or b(y), 
                xy = false;  
                x = 0; y = 0; %#ok<*NASGU> 
                % We need x and y for evaluation of string expressions  
                % like 'sin(x)+y'
                if b(2) == 0, 
                    lengthq = b(3);
                    lengthg = b(4);
                    % Neumann  BC on this edge segment       
                    try % checks if g = g(x,y) or g = const                 
                        % if it works, q and g are written in terms of (x,y)
                        a1 = eval(char(b(5:5+b(3)-1)));
                        a2 = eval(char(b(5+lengthq:5+lengthq+lengthg-1)));
                        if ~(isnumeric(a1)&&isnumeric(a2))
                            throw(finiteElement.wrongClass)
                        end
                    catch ME
                        % Tricky: If an exception was thrown, it must be
                        % parametrisied by x and y variables.
                        xy = true;                         
                    end
                else
                    % Dirichlet  BC on this edge segment
                    try
                        % if it works, h is written in terms of (x,y)
                        lengthh = b(5);
                        lengthr = b(6);
                        a1 = eval(char(b(9:9+lengthh-1)));
                        a2 = eval(char(b(9+lengthh:9+lengthh+lengthr-1)));
                        if ~(isnumeric(a1)&&isnumeric(a2))
                            throw(finiteElement.wrongClass)
                        end
                    catch ME 
                        xy = true; 
                    end
                end
              end        
        
            
            function val = isArclengthParametrized(b)
                % Checks if the function in b are arc-length-parametrized  
                % or not by trial-and-error.
                 
                val = false;  
                x = 0; y = 0; %#ok<*NASGU> 
                % We need x and y for evaluation of string expressions  
                % like 'sin(x)+y'
                if b(2) == 0
                    lengthq = b(3);
                    lengthg = b(4);
                    % Neumann  BC on this edge segment       
                    try % checks if g = g(x,y) or g = const                 
                        % if it works, q and g are written in terms of (x,y)
                        a1 = eval(char(b(5:5+b(3)-1)));
                        a2 = eval(char(b(5+lengthq:5+lengthq+lengthg-1)));
                        if ~(isnumeric(a1)&&isnumeric(a2))
                            throw(finiteElement.wrongClass)
                        end
                    catch % Don#t catch it
                        % Tricky: If an exception was thrown, it must be
                        % parametrisied in arclength, i.e. g = g(s) etc.
                        val = true;                         
                    end
                else
                    % Dirichlet  BC on this edge segment
                    try
                        % if it works, h is written in terms of (x,y)
                        lengthh = b(5);
                        lengthr = b(6);
                        a1 = eval(char(b(9:9+lengthh-1)));
                        a2 = eval(char(b(9+lengthh:9+lengthh+lengthr-1)));
                        if ~(isnumeric(a1)&&isnumeric(a2))
                            throw(finiteElement.wrongClass)
                        end
                    catch  % Don#t catch it 
                        % if it not works, it must be written as g(s)
                        val = true; 
                    end
                end
            end            
        end 
    end   
end
