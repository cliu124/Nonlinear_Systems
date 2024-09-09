classdef lagrange11D < finiteElements1D
    properties(Constant,Hidden)
        % stiff element
        S = reshape([1 -1;...
                    -1 1],4,1);
        % mass element
        M = reshape(1/6*[2 1;...
                         1 2],4,1);
        % convection
        C = reshape([-1/2 1/2;...
                     -1/2 1/2],4,1);
        % RHS element
        F = 0.5*[1 1]';    
        % index of neighbours
        % trivial but for generality.
        idx = 1:2;
        assembleStep = 1;
    end
    
    methods( Static,Access = public)  
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            idxvec0 = reshape(idx,1,np*2);
            idxvec1 = reshape([idx;idx],1,np*4);
            idxvec2 =  reshape([idx(1,:);idx(1,:);...
                                idx(2,:);idx(2,:)]...
                                                  ,1,4*np);     
        end 
        function [J] = makeJ(gridObj) 
            scl = superclasses(gridObj);
            for k = 1:length(scl)
                isgrid = strcmp(scl{k},'gridd');
                if isgrid
                    break;
                end
            end
            if ~isgrid
                throw(obj.wrongClass);
            end
            J = gridObj.p(2:end)-gridObj.p(1:end-1);           
        end
        
        function fluxThrougElementEdges= fluxThroughEdges(obj,u,c)
            % obj.fluxThroughEdges(grid,u,c)
            % Static method.
            % Arguments check:
            % 14, size(c), pause, %HU
            
            switch length(c) % c scalar! 
                case {1, obj.nElements}              
                    % c fits!
                case obj.nPoints
                    c = obj.point2Center(c);
                otherwise  % c is miss-dimensioned
                  %  c=ones(c(1,1),obj.nElements)*f; % HU                    
                   obj.wrongFormat.throwAsCaller;  
            end
            % 15, size(u), obj.nPoints % HU 
            if length(u)~=obj.nPoints
                %20,  u=[u; u(end)]; % HU 
               obj.wrongFormat.throwAsCaller; 
            end
            % elementArea is in 1D the same as dx!
            [~, elementArea] = obj.sideLengthAndArea;
            % 2 x nElements matrix, rather formal...
    
            fluxThrougElementEdges = [1;-1]*...
                (c.*(u(obj.t(2,:))-u(obj.t(1,:)))'./elementArea);
        end
               
        
        function jumps = fluxJumps(obj,fluxThroughElementEdges,~)
            % --- this is dimension depending -- here 1D sub intervals
            interiorEdges=sparse(obj.t([1 2],:),...
                obj.t([1 2],:),1,obj.nPoints,obj.nPoints);
            
            % intj+intj.' is 2 if interior edge and 1 if exterior edge
            interiorEdges=round(interiorEdges/3);
            jmps = sparse(obj.t([1 2],:),obj.t([1 2],:),...
                fluxThroughElementEdges([1 2],:),obj.nPoints,obj.nPoints); 
            jmps = interiorEdges.*abs(jmps);             
            jumps  = sum(jmps(obj.t(1,:),obj.t(2,:)).^2,1);  
        end
    end 
    methods(Access = public)
        function [DX,DY] = gradientMatrices(obj,grid)
            %%
            % * gradientMatrices IN: gridd OUT: double 
            % 
            % Method that computes matrix DX  , such that 
            %
            % grad u =  (DX*u)'                        
            %
            % at the center of all elements. 
            % DX is nElements x nPoints matrix. Note that this is
            % not an exact computation but an approximation with only
            % linear convergence. 
            
            % The indices of the first and second  points in each triangle. 
            idx1 = grid.t(1,:); idx2 = grid.t(2,:);                         
            % Compute Jacobi determinat  
            J = obj.makeJ(grid);              
            % FD 
            p21 = 1./J;           
            DX = sparse([1:grid.nElements,1:grid.nElements],...
                [idx1 idx2],[-p21 p21] ,...
                grid.nElements,grid.nPoints);            
             
        end
    end
end