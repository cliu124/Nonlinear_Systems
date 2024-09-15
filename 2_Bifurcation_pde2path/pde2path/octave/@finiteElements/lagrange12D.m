classdef lagrange12D < finiteElements2D
    % class for build Lagrange L1 FE in 2D
    % Inherits everything from abstract class finiteElements
    % Only different Matrices were defined as constant props.
        
    properties(Constant)
        % local matrices, resp. vectors. We use the vector form to use
        % later the sparse function to create the matrices... cf. the code
        % in assema etc.
        
        % stiff
        S1 = 0.5*[ 1
                  -1
                   0
                  -1
                   1
                   0
                   0
                   0
                   0];
               
        S2 = 0.5*[ 2 
                  -1 
	              -1 
                  -1 
                   0 
                   1  
                  -1
                   1
                   0];
              
        S3 = 0.5*[ 1  
                   0  
	              -1  
                   0 
                   0  
                   0  
                  -1
                   0
                   1];
              
        % mass
        M = 1/24*[2  
	               1  
	               1  
                   1   
                   2  
                   1  
                   1
                   1 
                   2];
               
        % RHS vector
        F = 1/6*[1
                  1
                  1]; 
              
        % Matrices convection terms 
        % for use by "sparse" in vector form
        C1 = [-1/6   
           -1/6   
           -1/6  
            1/6
            1/6
            1/6
            0
            0
            0];
        
        C2 = [-1/6   
           -1/6   
           -1/6 
            0
            0
            0
            1/6
            1/6
            1/6]; 
        
      
        % sparsity constant for allocating memory
        sparsityConstant = 5;
        
        % index vector for selecting coordinates
        idx = 1:3;
    end
    
    methods(Static,Access=public)   
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            idxvec0 = reshape(idx,1,np*3);
            idxvec1 = reshape([idx;idx;idx],1,np*9);
            idxvec2 =  reshape([idx(1,:);idx(1,:);idx(1,:);...
                idx(2,:); idx(2,:); idx(2,:);...
                idx(3,:);idx(3,:);idx(3,:)],1,9*np);     
        end  
     %Implement abstract inheritance from finiteElements class
        function ddncu = fluxThroughEdges(obj,u,c)
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
            for k = 1:3,
                k1 = rem(k ,3)+1;
                k2 = rem(k1,3)+1;
                dx(k,:) = obj.p(1,obj.t(k1,:)) - obj.p(1,obj.t(k2,:));
                dy(k,:) = obj.p(2,obj.t(k1,:)) - obj.p(2,obj.t(k2,:));
            end;
            % gradients of basis functions
            g1x=0.5*dy(1,:)./elementArea;  g2x=0.5*dy(2,:)./elementArea;
            g3x=0.5*dy(3,:)./elementArea;  g1y=-0.5*dx(1,:)./elementArea;
            g2y=-0.5*dx(2,:)./elementArea; g3y=-0.5*dx(3,:)./elementArea;

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
    methods(Access=public)
    function [DX,DY] = gradientMatrices(obj,grid)
        %%
        % * gradientMatrices IN: gridd OUT: double,double
        % Method that computes matrices DX DY, such that 
        %
        % grad u = [(DX*u)' ,(DY*u)']                       
        %
        % at the center of all triangles. 
        % DX and DY are nElements x nPoints matrices. 
        
        % We want to use makeJ, so it is not static...
        
        % The indices of the first, second and third points in each
        % triangle. 
        idx1 = grid.t(1,:); idx2 = grid.t(2,:); idx3 = grid.t(3,:);
        
        % p1 are the  values of the "first" point, p2 the values of
        % the "second" point p3 the values of the third point 
        % of all triangles. (2 x nElement matrices)
        p1 = grid.p(:,(idx1)); p2 = grid.p(:,(idx2)); p3 = grid.p(:,(idx3));
        
        % Compute Jacobi determinat  
        J = obj.makeJ(grid);
        
        p12 = (p1(2,:)-p2(2,:))./J; p31 = (p3(2,:)-p1(2,:))./J;
        p23 = (p2(2,:)-p3(2,:))./J;
        
        DX = sparse([1:grid.nElements,1:grid.nElements,1:grid.nElements],...
            [idx1 idx2 idx3],[p23 p31 p12] ,...
            grid.nElements,grid.nPoints);
        
        p21 = (p2(1,:)-p1(1,:))./J;  p13 = (p1(1,:)-p3(1,:))./J;
        p32 = (p3(1,:)-p2(1,:))./J;
        
        DY = sparse([1:grid.nElements,1:grid.nElements,1:grid.nElements],...
            [idx1 idx2 idx3],[p32 p13 p21] ,...
            grid.nElements,grid.nPoints);
    end
    end
end