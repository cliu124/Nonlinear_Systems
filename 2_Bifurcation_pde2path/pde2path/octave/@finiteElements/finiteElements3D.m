classdef  finiteElements3D < finiteElements
    % Abstract FiniteElements 3D class definition
    % Note: assema is defined in Superclass finiteElements
    
    methods(Access = public)
        function [K,M,F] = createMatrixEntries(obj,gridObj,cf,af,ff)
            % the sizes of elementary matrices, varies from element
            % type to element type...
            
            sizeVector = sqrt(size(obj.S1,1));
            sizeMatrix = size(obj.S1,1);                      
            
            % compute the coefficients. The result is always a k x nElelemts matrix
            % where k = 3 or 9 for c. k is 3 if and only if the
            % diffussioncoefficient
            % is a diagonal matrix or scalar.
            % k = 1 for a and f.
            [cval,aval,fval] = obj.aCoefficients(gridObj,cf,af,ff);
            
            [K1,K2,K3,K4,K5,K6,K7,K8,K9] = obj.getConstantPartOfStiffnessMatrix(gridObj);                       
           switch size(cval,1) 
                case 1
                    K = reshape(K1*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(1,:))...
                    + K5*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(1,:))...
                    + K9*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(1,:))...
                        ,1,gridObj.nElements*sizeMatrix);
                case 3
                K = reshape(K1*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(1,:))...
                    + K5*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(2,:))...
                    + K9*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(3,:))...
                        ,1,gridObj.nElements*sizeMatrix);
                case 9 
                K = reshape(K1*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(1,:))...
                    + K2*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(2,:))...
                    + K3*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(3,:))...
                    + K4*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(4,:))...
                    + K5*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(5,:))...
                    + K6*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(6,:))...
                    + K7*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(7,:))...
                    + K8*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(8,:))...
                    + K9*sparse(1:gridObj.nElements,1:gridObj.nElements,cval(9,:))...
                        ,1,gridObj.nElements*sizeMatrix);
                otherwise
                    error(['The object cval has ' num2str(size(cval,1)) ' Rows.'...
                        ' If Diffusion Coefficient is scalar, cval must have three equal rows. '...
                        'If it is a diagonal matrix, cval must have these three diagonalelements as rows. '...
                        'If it is a full matrix, cval must have nine rows, which are taken from the matrix columnwise. '...
                        'So if the rows are labeld with c1,c2,c3,c4, and so on, then the matrix entries are'...
                        ' [c1 c4 c7; c2 c5 c8; c3 c6 c9].'...
                        ' Other syntaxes are not supported. For further help, see aCoefficients() in your used grid class.']);
            end
            
            M = reshape(obj.M*(aval.*obj.makeJ(gridObj)),1,gridObj.nElements*sizeMatrix);
            
            F = reshape(obj.F*(obj.makeJ(gridObj).*fval),1,gridObj.nElements*sizeVector);
        end
        
        function val = createConvectionEntries(obj,gridObj,b)
           
            sizeMatrix = size(obj.S1,1);
           
            
            p1 = gridObj.p(:,(gridObj.t(1,:)));
            p2 = gridObj.p(:,(gridObj.t(2,:)));
            p3 = gridObj.p(:,(gridObj.t(3,:)));
            p4 = gridObj.p(:,(gridObj.t(4,:)));
            
            % get element corner coordinates
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:); 
            x41 = p4(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:); 
            y41 = p4(2,:)-p1(2,:); 
            z21 = p2(3,:)-p1(3,:);
            z31 = p3(3,:)-p1(3,:); 
            z41 = p4(3,:)-p1(3,:);
            
            % now we compute J for every element in one operation 
            
            % Here we must use the "analytic" formulas...
            
             
            
            xi_x = (y31.*z41-y41.*z31) ;
            eta_x = (y41.*z21-y21.*z41) ;
            zeta_x =  (y21.*z31-y31.*z21) ;
            
            xi_y = (x41.*z31-x31.*z41) ;
            eta_y = (x21.*z41-x41.*z21);
            zeta_y = (x31.*z21-x21.*z31);
            
            xi_z = (x31.*y41-x31.*y41);
            eta_z = (x41.*y21-x21.*y41);
            zeta_z = (x21.*y31-x31.*y21);
            
            
            % to be implemented
            bval = obj.convCoefficients(gridObj,b);
         
            J1 = (bval(1,:).*xi_x+bval(2,:).*xi_y+bval(3,:).*xi_z);
            J2 = (bval(1,:).*eta_x+bval(2,:).*eta_y+bval(3,:).*eta_z);
            J3 = (bval(1,:).*zeta_x+bval(2,:).*zeta_y+bval(3,:).*zeta_z);
             
            val = reshape(obj.C1*J1+obj.C2*J2+obj.C3*J3,1,gridObj.nElements*sizeMatrix);            
        end
   

    end
    
    methods(Static,Access=public)
        function [cval,aval,fval] =  aCoefficients(gridObj,cc,aa,ff)
            [cval,aval,fval] = gridObj.aCoefficientsMpt(cc,aa,ff);
        end 
        
        function cval = convCoefficients(gridObj,b)
            cval = gridObj.cCoefficients(b);
        end
          
        function J = makeJ(gridObj)
                % helper method to compute the Jacobian determinant
                %                        
                
                if ~(isa(gridObj,'grid3D')||isa(gridObj,'grid3Dpr'))
                    gridObj.wrongClass.throwAsCaller;
                end
                p1 = gridObj.p(:,(gridObj.t(1,:)));
                p2 = gridObj.p(:,(gridObj.t(2,:)));
                p3 = gridObj.p(:,(gridObj.t(3,:)));
                p4 = gridObj.p(:,(gridObj.t(4,:)));
                % we don't need p5-p6

                % get element corner coordinates
                x21 = p2(1,:)-p1(1,:);
                x31 = p3(1,:)-p1(1,:); 
                x41 = p4(1,:)-p1(1,:);
                y21 = p2(2,:)-p1(2,:);
                y31 = p3(2,:)-p1(2,:); 
                y41 = p4(2,:)-p1(2,:); 
                z21 = p2(3,:)-p1(3,:); % zero for prisms
                z31 = p3(3,:)-p1(3,:); % zero for prisms
                z41 = p4(3,:)-p1(3,:);

                % now we compute J for every element in one operation 

                % Here we must use the "analytic" formulas...
                % The Determinat of Jacobian
                J =  (x21.*y31-x31.*y21).*z41...
                        -(x21.*y41-x41.*y21).*z31...
                        +(x31.*y41-x41.*y31).*z21;               

        end
    end
    
    methods(Access = public, Hidden)               
        function [K11,K21,K31,K12,K22,K32,K13,K23,K33] = getConstantPartOfStiffnessMatrix(obj,gridObj)
            % helper method to compute the former abstract Js
            %                        
            p1 = gridObj.p(:,(gridObj.t(1,:)));
            p2 = gridObj.p(:,(gridObj.t(2,:)));
            p3 = gridObj.p(:,(gridObj.t(3,:)));
            p4 = gridObj.p(:,(gridObj.t(4,:)));
            % we don't need p5-p6
            
            % get element corner coordinates
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:); 
            x41 = p4(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:); 
            y41 = p4(2,:)-p1(2,:); 
            z21 = p2(3,:)-p1(3,:); % zero for prisms
            z31 = p3(3,:)-p1(3,:); % zero for prisms
            z41 = p4(3,:)-p1(3,:);
            
            % now we compute J for every element in one operation 
            
            % Here we must use the "analytic" formulas...
            % The Determinat of Jacobian
            J =  ((x21.*y31-x31.*y21).*z41...
                    -(x21.*y41-x41.*y21).*z31...
                    +(x31.*y41-x41.*y31).*z21);
           
           
            
            % 
            % Inversion by using the famous 3x3 inversion formula
            xi_x = (y31.*z41-y41.*z31); 
            eta_x = (y41.*z21-y21.*z41);
            zeta_x =  (y21.*z31-y31.*z21); % zero for prisms
            
            xi_y = (x41.*z31-x31.*z41) ;
            eta_y = (x21.*z41-x41.*z21);  
            zeta_y = (x31.*z21-x21.*z31); % zero for prisms
            
            xi_z = (x31.*y41-x41.*y31);
            eta_z = (x41.*y21-x21.*y41);
            zeta_z = (x21.*y31-x31.*y21);
            
            K11 = obj.S1*((xi_x.^2)./J) + obj.S2*((eta_x.^2)./J)...
                + obj.S3*((zeta_x.^2)./J) + (obj.S7 +obj.S4)*((xi_x.*eta_x)./J)...
                + (obj.S8 + obj.S5)*((xi_x.*zeta_x)./J) + (obj.S9 + obj.S6)*((eta_x.*zeta_x)./J);
            
            K22 = obj.S1*((xi_y.^2)./J) + obj.S2*((eta_y.^2)./J)...
                + obj.S3*((zeta_y.^2)./J) + (obj.S7 +obj.S4)*((xi_y.*eta_y)./J)...
                + (obj.S8 + obj.S5)*((xi_y.*zeta_y)./J) + (obj.S9 + obj.S6)*((eta_y.*zeta_y)./J);
            
            K33 = obj.S1*((xi_z.^2)./J) + obj.S2*((eta_z.^2)./J)...
                + obj.S3*((zeta_z.^2)./J) + (obj.S7 +obj.S4)*((xi_z.*eta_z)./J)...
                + (obj.S8 + obj.S5)*((xi_z.*zeta_z)./J) + (obj.S9 + obj.S6)*((eta_z.*zeta_z)./J);  
            
            K21 = obj.S1*((xi_x.*xi_y)./J) + obj.S2*((eta_x.*eta_y)./J)...
                + obj.S3*((zeta_x.*zeta_y)./J) + obj.S4*((xi_y.*eta_x)./J)...
                + obj.S5*((xi_y.*zeta_x)./J) + obj.S6*((eta_y.*zeta_x)./J)...
                + obj.S7*((xi_x.*eta_y)./J)...
                + obj.S8*((xi_x.*zeta_y)./J) + obj.S9*((eta_x.*zeta_y)./J);
                        
            K12 = obj.S1*((xi_x.*xi_y)./J) + obj.S2*((eta_x.*eta_y)./J)...
                + obj.S3*((zeta_x.*zeta_y)./J) + obj.S4*((xi_x.*eta_y)./J)...
                + obj.S5*((xi_x.*zeta_y)./J) + obj.S6*((eta_x.*zeta_y)./J)...
                + obj.S7*((xi_y.*eta_x)./J)...
                + obj.S8*((xi_y.*zeta_x)./J) + obj.S9*((eta_y.*zeta_x)./J);
                
            K31 = obj.S1*((xi_x.*xi_z)./J) + obj.S2*((eta_x.*eta_z)./J)...
                + obj.S3*((zeta_x.*zeta_z)./J) + obj.S4*((xi_z.*eta_x)./J)...
                + obj.S5*((xi_z.*zeta_x)./J) + obj.S6*((eta_z.*zeta_x)./J)...
                + obj.S7*((xi_x.*eta_z)./J)...
                + obj.S8*((xi_x.*zeta_z)./J) + obj.S9*((eta_x.*zeta_z)./J);
            
            K13 = obj.S1*((xi_x.*xi_z)./J) + obj.S2*((eta_x.*eta_z)./J)...
                + obj.S3*((zeta_x.*zeta_z)./J) + obj.S4*((xi_x.*eta_z)./J)...
                + obj.S5*((xi_x.*zeta_z)./J) + obj.S6*((eta_x.*zeta_z)./J)...
                + obj.S7*((xi_z.*eta_x)./J)...
                + obj.S8*((xi_z.*zeta_x)./J) + obj.S9*((eta_z.*zeta_x)./J);
                       
            K32 = obj.S1*((xi_z.*xi_y)./J) + obj.S2*((eta_z.*eta_y)./J)...
                + obj.S3*((zeta_z.*zeta_y)./J) + obj.S4*((xi_z.*eta_y)./J)...
                + obj.S5*((xi_z.*zeta_y)./J) + obj.S6*((eta_z.*zeta_y)./J)...
                + obj.S7*((xi_y.*eta_z)./J)...
                + obj.S8*((xi_y.*zeta_z)./J) + obj.S9*((eta_y.*zeta_z)./J);
            
            K23 = obj.S1*((xi_z.*xi_y)./J) + obj.S2*((eta_z.*eta_y)./J)...
                + obj.S3*((zeta_z.*zeta_y)./J)+ obj.S4*((xi_y.*eta_z)./J)...
                + obj.S5*((xi_y.*zeta_z)./J) + obj.S6*((eta_y.*zeta_z)./J)...
                 + obj.S7*((xi_z.*eta_y)./J)...
                + obj.S8*((xi_z.*zeta_y)./J) + obj.S9*((eta_z.*zeta_y)./J);   
        end
       
    end  
        
    methods(Static,Access=public)       
        % [Q,G,H,R] = assemb(gridObj)  
        % To be overwritten in derived classes        
    end
    
    methods(Static,Access = public)
        % fluxThrougElementEdges= fluxThroughEdges(obj,u,c)  
        % normFminusau = localErrorL2(obj,a,f)  
        % jumps = fluxJumps(obj,fluxThroughElementEdges,order)
    end  
end

