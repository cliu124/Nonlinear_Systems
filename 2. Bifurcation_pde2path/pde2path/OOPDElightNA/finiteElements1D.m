classdef finiteElements1D < finiteElements
    methods(Static, Access = public) 
        function [Q,G,H,R] = assemb(gridObj)
            % to fulfill the signature of all assemb functions,
            % it must be non static, despite we do not need obj here...
            b = gridObj.b; 
            rE =  gridObj.e(1,2);
            n = gridObj.nPoints;
            
            Q = sparse(n, n); H =  sparse(n, n); R = sparse(n,1);G = sparse(n,1);
            
            % left BD
            m = b(2,1); lengthq = b(3,1); lengthg = b(4,1);
            
            if m == 0 % case Robin
                Q(1,1) = eval(char(b(5:5+lengthq-1,1)));
                G(1) = eval(char(b(5+lengthq:5+lengthq+lengthg-1,1)));
            else % case Dirichlet BC
                lengthh = b(5,1); lengthr = b(6,1);
                Q(1,1) = eval(char(b(7:7+lengthq-1,1)));
                G(1) = eval(char(b(7+lengthq:7+lengthq+lengthg-1,1)));                
                H(1,1) = eval(char(b(9:9+lengthh-1,1)));
                R(1) = eval(char(b(9+lengthh:9+lengthh+lengthr-1,1)));
            end             
            % right BD            
            m = b(2,2); lengthq = b(3,2); lengthg = b(4,2);
            
            if m == 0 % case Robin
                Q(rE,rE) = eval(char(b(5:5+lengthq-1,2)));
                G(rE) = eval(char(b(5+lengthq:5+lengthq+lengthg-1,2)));
                H(rE,rE)=0;
                R(rE)=0;
            else % case   Dirichlet BC
                lengthh = b(5,2);
                lengthr = b(6,2);
                Q(rE,rE) = eval(char(b(7:7+lengthq-1,2)));
                G(rE) = eval(char(b(7+lengthq:7+lengthq+lengthg-1,2)));
                H(rE,rE) = eval(char(b(9:9+lengthh-1,2)));
                R(rE) = eval(char(b(9+lengthh:9+lengthh+lengthr-1,2)));
            end
        end         
    end  
    
    methods(Access = public)
        function [K,M,F] = createMatrixEntries(obj,gridObj,cf,af,ff) 
            sizeVector = sqrt(size(obj.S,1));
            sizeMatrix = size(obj.S,1);
            nt = gridObj.nElements;
            [cval,aval,fval] = obj.aCoefficients(gridObj,cf,af,ff);
            J = obj.makeJ(gridObj);  
            
            K1 = getConstantPartOfStiffnessMatrix(obj,gridObj);
           
            K = reshape(K1*sparse(1:gridObj.nElements,...
                1:gridObj.nElements,cval(1,:))...
                     ,1,gridObj.nElements*sizeMatrix);
        
            M = reshape(obj.M*(aval.*J),1,nt*sizeMatrix);            
            F = reshape(obj.F*(J.*fval),1,nt*sizeVector);            
        end
        
        function val = createConvectionEntries(obj,gridObj,b)  
            sizeMatrix = size(obj.C,1);
            J = obj.convCoefficients(gridObj,b); 
            val = reshape(obj.C*J,1,gridObj.nElements*sizeMatrix);        
        end       
    end
    methods(Access=public,Hidden)
        function K1 = getConstantPartOfStiffnessMatrix(obj,gridObj)
            J = obj.makeJ(gridObj);             
            K1 = obj.S*(1./J);
        end
    end
    methods(Static,Access=protected)       
        
        function [cval,aval,fval] =  aCoefficients(gridObj,cc,aa,ff)
            [cval,aval,fval] = gridObj.aCoefficientsMpt(cc,aa,ff);
        end
        
        function bval = convCoefficients(gridObj,b)
            bval = gridObj.convCoefficientsMpt(b);
        end
    end
    
    methods(Static,  Access = public)
      [J] = makeJ(gridObj) 
    end
end

