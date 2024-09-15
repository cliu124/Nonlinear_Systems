classdef lagrange02D < finiteElements2D
    %class stub that implements L0 elements
     
    
    properties(Constant)       
        idx = [1,2];
    end
     
    
    methods(Static,Access = public)
        function M = mass(grid)
            j =  lagrange02D.makeJ(grid);            
            M = 0.5*sparse(1:grid.nElements,...
                1:grid.nElements,j);  
        end
        function makeIndex()
        end
    end
    
end

