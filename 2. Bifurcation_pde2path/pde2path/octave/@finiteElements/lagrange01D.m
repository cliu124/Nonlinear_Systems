classdef lagrange01D < finiteElements1D
    %class stub that implements L0 elements
     
    
    properties(Constant)
        
        idx = [1,2];
    end
     
    
    methods(Static,Access = public)
        function M = mass(grid)
            j = grid.p(1,2:end)-grid.p(1,1:end-1);            
            M = sparse(1:grid.nElements,...
                1:grid.nElements,j);  
        end
        function makeIndex()
        end
    end
    
end