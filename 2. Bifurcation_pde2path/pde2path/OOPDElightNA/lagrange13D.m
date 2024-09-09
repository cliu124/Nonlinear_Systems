classdef lagrange13D< finiteElements3D
    % lagrange13D class definition
    
    properties(Constant)
        %  stiff
        S1 =  reshape([1/6    -1/6     0     0
                      -1/6     1/6     0     0  
                       0     0     0     0
                       0     0     0     0] ,16,1);        
        
        S2 = reshape( [1/6     0    -1/6     0
                       0     0     0     0
                      -1/6     0     1/6     0
                       0     0     0     0],16,1);

        S3 = reshape( [1/6     0     0    -1/6
                       0     0     0     0
                       0     0     0     0
                      -1/6     0     0     1/6],16,1);
                  
                  
        % Note: Here we have to multiply the integral by the
        % factor of two, if we want so include the
        % symmetrie of intefgral here and not in the
        % assembling... cf.  method createMatrixEntries of
        % finiteElements3D class
        % But the integral is not symmetric as long as c is not symmetric.
        % So we do not multiply by 2.
      

        S4 =  reshape( [ 1     0    -1     0
                        -1     0     1     0
                         0     0     0     0
                         0     0     0     0 ] /6 ,16,1);
       
       
        S5 =  reshape( [ 1     0     0    -1
                        -1     0     0     1
                         0     0     0     0
                         0     0     0     0 ] / 6 ,16,1  );


        S6 =  reshape( [1     0     0    -1
                         0     0     0     0
                        -1     0     0     1
                         0     0     0     0] / 6,16,1);

                  
        % S7 to S9 represents the non symetric case. So we do not have to
        % multiply the integrals by 2 because we do not regard c as
        % symmetric
        
        
        S7 =  reshape( [1    -1     0     0
                         0     0     0     0
                        -1     1     0     0
                         0     0     0     0 ] / 6 ,16,1);
                     
        S8 =  reshape( [ 1    -1     0     0
                         0     0     0     0
                         0     0     0     0
                        -1     1     0     0 ] / 6 ,16,1  );


        S9 =  reshape( [1     0    -1     0
                         0     0     0     0
                         0     0     0     0
                        -1     0     1     0] / 6  ,16,1);

                     
        
     

 
                  
 
       
       % Element Mass      
       
        M = reshape([1 / 60      1 / 120	 1 / 120	 1 / 120	
                     1 / 120	 1 / 60	     1 / 120	 1 / 120	
                     1 / 120	 1 / 120	 1 / 60	     1 / 120	
                     1 / 120	 1 / 120	 1 / 120	 1 / 60],16,1);  
 
        % convection 3x
        C1= reshape(1/24*[  -1 1 0 0
                            -1 1 0 0
                            -1 1 0 0
                            -1 1 0 0],16,1);
        
        
        C2 =  reshape(1/24*[ -1 0 1 0
                             -1 0 1 0
                             -1 0 1 0
                             -1 0 1 0],16,1);
 
         
        C3 =  reshape(1/24*[ -1 0 0 1
                             -1 0 0 1
                             -1 0 0 1
                             -1 0 0 1],16,1);
          
        % RHS vector  
        F = 1/24 * [1
                    1
                    1
                    1];
           
           
        % the points numbers, we have linear functions on
        % tetrahedrals, so we need only four points...
        idx = 1:4;
    end 
    
    methods(Static,Access = public)   
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            % reshape index vectors for using sparse in assema
            idxvec0 = reshape(idx,1,np*4);
            idxvec1 = reshape([idx;idx;idx;idx],1,np*16);
            idxvec2 =  reshape([idx(1,:);idx(1,:);idx(1,:);idx(1,:);...
                idx(2,:); idx(2,:); idx(2,:);idx(2,:);...
                idx(3,:);idx(3,:);idx(3,:);idx(3,:);...
                idx(4,:);idx(4,:);idx(4,:);idx(4,:)],1,16*np);        
        end    
    end
    
    methods(Static,Access=public)
        function [Q,G,H,R] = assemb(grid)
            % [Q,G,H,R] = assemb(grid)
            % Static method for   3D terahedral elements
            % assembles the boundary condition matrix.
            % 
            % (c) 2013 by Uwe Prüfert
            %
            % ChangeLog
            % 03/03/2015 by Tobias Höer
            % fixed a bug by changing Line from
            %'p22 = [a;zeros(1,grid.nEdges)];' to
            %'p22 = [b;zeros(1,grid.nEdges)];'



            nPts = grid.nPoints;
            nEdges = grid.nEdges;
             
            qval = zeros(1,grid.nEdges);
            gval = zeros(1,grid.nEdges);
            hval = zeros(1,grid.nEdges);
            rval = zeros(1,grid.nEdges);
 
            for k = 1:grid.nEdges   
                    [qval(k),gval(k),hval(k),rval(k)] = grid.bCoefficients(...
                        grid.p(:,grid.e(1:3,k)),grid.b(:,grid.e(5,k)));
            end
        
             
           
            %  We have to transform the triangles from the boundary
            %  from 3D into 2D. We use elementary calculations. Idea:
            %  Set p11 = (0,0) and p22 = (||p2-p1||,0). We calculate
            %  p33 by using \cos \alpha = \frac{b^2 + c^2 - a^2}{2bc}

            p1 = grid.p(1:3,grid.e(1,:));
            p2 = grid.p(1:3,grid.e(2,:));
            p3 = grid.p(1:3,grid.e(3,:));
            
            a = sqrt(sum((p2-p1).^2,1));
            b = sqrt(sum((p3-p2).^2,1));
            c = sqrt(sum((p3-p1).^2,1));
            
            
            
            alpha = acos((b.^2+c.^2-a.^2)./2./b./c);
            
            p11 = [0;0];
            p22 = [b;zeros(1,grid.nEdges)];
            p33 = [cos(alpha).*c;sin(alpha).*c];
            
            % get element corner coordinates
            x21 = p22(1,:)-p11(1,:);
            x31 = p33(1,:)-p11(1,:);
            y21 = p22(2,:)-p11(2,:);    
            y31 = p33(2,:)-p11(2,:);

            J = (x21.*y31-x31.*y21);
            
            

           
            % in 3D tetrehedral elment case, the boundary segments are
            % triangles. We can use the 2D Mass and RHS matrix/vector.
            Qe = reshape(lagrange12D.M*(J.*qval),1,9*nEdges);
            ge = reshape(lagrange12D.F*(J.*gval),1,3*nEdges); 
 
            He = reshape([1 0 0 0 1 0 0 0 1]'*(J.*hval),1,9*nEdges);  
            re = reshape([1 1 1]'*(J.*rval),1,3*nEdges); 

            % re-aranging like in assema...
            indxTriPts = grid.e(1:3,:);
            % now we use the 2D stuff to compute the index.structure
            [indx0,indx1,indx2] = lagrange12D.makeIndex(indxTriPts,nEdges);        

            Q = sparse(indx1,indx2,Qe,nPts,nPts); 
            G = sparse(indx0,1,ge,nPts,1);           
            H = sparse(indx1,indx2,He,nPts,nPts);
            
             
            R = sparse(indx0,1,re,nPts,1);          
        end 
    end
    methods(Access=public)
function [DX,DY,DZ] = gradientMatrices(obj,grid)
    %%
    % * gradientMatrices IN: gridd OUT: double,double
    % 
    % Method that computes matrices DX DY, such that 
    %
    % grad u = [(DX*u)' ,(DY*u)' (DZ*u)']                       
    %
    % at the center of all triangles. 
    % DX and DY are nElements x nPoints matrices. Note that this is
    % not an exact computation but an approximation with only
    % linear convergence. 
    
    % We want to use makeJ, so it is not static...
    
    % The indices of the first, second and third points in each
    % triangle. 
    idx1 = grid.t(1,:); idx2 = grid.t(2,:); idx3 = grid.t(3,:);
    idx4 = grid.t(4,:);
    
    % p1 are the  values of the "first" point, p2 the values of
    % the "second" point, etc. (2 x nElement matrices)
    p1 = grid.p(:,(idx1)); p2 = grid.p(:,(idx2));
    p3 = grid.p(:,(idx3)); p4 = grid.p(:,(idx4));
    
    % Compute Jacobi determinat (is it correct?)
    J = obj.makeJ(grid);
    
    % x_2-x_1 etc. in every element
    x21 = (p2(1,:)-p1(1,:)); x31 = (p3(1,:)-p1(1,:)); x41 = (p4(1,:)-p1(1,:));
    
    % y_2-y_1 etc. in every element
    y21 = (p2(2,:)-p1(2,:)); y31 = (p3(2,:)-p1(2,:)); y41 = (p4(2,:)-p1(2,:));
        
    % z_2-z_1 etc. in every element
    z21 = (p2(3,:)-p1(3,:)); z31 = (p3(3,:)-p1(3,:)); z41 = (p4(3,:)-p1(3,:));
    
    %     | x21 y21 z21 |
    % A = | x31 y31 z31 |
    %     | x41 y41 z41 |
    %
    %  -1       |y31*z41-z31*y41  z21*y41-y21*z41 y21*z31-z21*y31 |
    % A  =  1   |z31*x41-x31*z41  x21*z41+z21*x41 z21*x31-x21*z31 |
    %     det A |x31*y41-y31*x41  y21*x41+x21*y41 x21*y31-y21*x31 |
    %  
    %      | A B C |  
    %    = | D E F |
    %      | G H I |
    
    
    A = (y31.*z41-z31.*y41)./J;    
    B = (z21.*y41-y21.*z41)./J;
    C = (y21.*z31-z21.*y31)./J;
    D = (z31.*x41-x31.*z41)./J;
    E = (x21.*z41-z21.*x41)./J;
    F = (z21.*x31-x21.*z31)./J; %#ok<*PROPLC>
    G = (x31.*y41-y31.*x41)./J;
    H = (y21.*x41-x21.*y41)./J;
    I = (x21.*y31-y21.*x31)./J;
    
    DX = sparse([1:grid.nElements,1:grid.nElements,...
        1:grid.nElements,1:grid.nElements],...
        [idx1 idx2 idx3 idx4],[-(A+B+C) A B C] ,...
        grid.nElements,grid.nPoints);
    
    
    DY = sparse([1:grid.nElements,1:grid.nElements,...
        1:grid.nElements,1:grid.nElements],...
        [idx1 idx2 idx3 idx4],[-(D+E+F) D E F] ,...
        grid.nElements,grid.nPoints);
    
    DZ = sparse([1:grid.nElements,1:grid.nElements,...
        1:grid.nElements,1:grid.nElements],...
        [idx1 idx2 idx3 idx4],[-(G+H+I) G H I] ,...
        grid.nElements,grid.nPoints);
    end
    end
end

