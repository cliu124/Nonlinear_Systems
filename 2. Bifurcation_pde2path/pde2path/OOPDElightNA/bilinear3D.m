classdef bilinear3D < finiteElements3D
    % class definition for bilinear FEs on prism elements
    % in 3D
    % fem = bilinear3D();
    % (c) 2013, Uwe Prüfert
    
    % Hints: 
    %   * Ensure in assema etc. that the grid is a grid3Dpr  object.
    %   * The local matrices are in vector form. To make the code more
    %     readable, we write down the matrix-form and reshape it. This is a
    %     preparation for the use of the sparse function when assemble
    %     the matrices.
    %   * The indx must be indx = 1:6

    properties(Constant) 
        % local stiffness matrices
        
        % They are ordered in the following way:
        %{ S1 S4 S6;
        %  S7 S2 S5;
        %  S9 S8 S3 }
        
        S1 = reshape([ 1/6      -1/6	 0       1/12	-1/12	 0	
                      -1/6       1/6	 0      -1/12	 1/12	 0	
                        0        0       0       0       0       0	
                        1/12    -1/12	 0       1/6	-1/6	 0	
                       -1/12     1/12	 0      -1/6	 1/6	 0	
                        0        0       0       0	     0       0],...
                        36,1);
                   
        S2 = reshape([  1/6      0	-1/6	 1/12	 0	-1/12	
                        0        0	 0       0       0	 0	
                       -1/6	 0	 1/6	-1/12	 0	 1/12	
                        1/12	 0	-1/12	 1/6	 0	-1/6	
                        0        0	 0       0       0	 0	
                       -1/12	 0	 1/12	-1/6	 0	 1/6 ],...
                        36,1);
                  
        S3 = reshape([ 1/12	 1/24	 1/24	-1/12	-1/24	-1/24	
                       1/24	 1/12	 1/24	-1/24	-1/12	-1/24	
                       1/24	 1/24	 1/12	-1/24	-1/24	-1/12	
                      -1/12	-1/24	-1/24	 1/12	 1/24	 1/24	
                      -1/24	-1/12	-1/24	 1/24	 1/12	 1/24	
                      -1/24	-1/24	-1/12	 1/24	 1/24	 1/12 ],...
                      36,1);
               
                 
        S4 = reshape([       4     0    -4     2     0    -2
                            -4     0     4    -2     0     2
                             0     0     0     0     0     0
                             2     0    -2     4     0    -4
                            -2     0     2    -4     0     4
                             0     0     0     0     0     0  ] /24,...
                      36,1);
                  
         
        S5 = reshape([ 2     2     2    -2    -2    -2
                        -2    -2    -2     2     2     2
                         0     0     0     0     0     0
                         2     2     2    -2    -2    -2
                        -2    -2    -2     2     2     2
                         0     0     0     0     0     0]/24,...
                                          36,1); 
        S6 = reshape([   2     2     2    -2    -2    -2
                         0     0     0     0     0     0
                        -2    -2    -2     2     2     2
                         2     2     2    -2    -2    -2
                         0     0     0     0     0     0
                        -2    -2    -2     2     2     2]/24,...
                      36,1);
                 
        S7 = reshape([       4    -4     0     2    -2     0
                             0     0     0     0     0     0
                            -4     4     0    -2     2     0
                             2    -2     0     4    -4     0
                             0     0     0     0     0     0
                            -2     2     0    -4     4     0  ] /24,...
                      36,1);
                  
         
        S8 = reshape([   2    -2     0     2    -2     0
                         2    -2     0     2    -2     0
                         2    -2     0     2    -2     0
                        -2     2     0    -2     2     0
                        -2     2     0    -2     2     0
                        -2     2     0    -2     2     0]/24,...
                                          36,1); 
                                      
        S9 = reshape([   2     0    -2     2     0    -2
                         2     0    -2     2     0    -2
                         2     0    -2     2     0    -2
                        -2     0     2    -2     0     2
                        -2     0     2    -2     0     2
                        -2     0     2    -2     0     2]/24,...
                      36,1);
        
        
        %   convection
        
        % The old matrixes. They seem to be wrong.
        
%         C1 = reshape([  -1/18	0       -1/36	-1/36	 0      -1/72	
%                          0      1/18	 1/36	 0       1/36	 1/72	
%                         -1/36	1/36	 0      -1/72	 1/72	 0	
%                         -1/36	0       -1/72	-1/18	 0      -1/36	
%                          0      1/36	 1/72	 0       1/18	 1/36	
%                         -1/72	1/72	 0      -1/36	 1/36	 0  ],...
%                         36,1);
%                   
%          
%         C2  = reshape([  -1/18	-1/36	 0      -1/36	-1/72	 0	
%                         -1/36	 0       1/36	-1/72	 0       1/72	
%                          0       1/36	 1/18	 0       1/72	 1/36	
%                         -1/36	-1/72	 0      -1/18	-1/36	 0	
%                         -1/72	 0       1/72	-1/36	 0       1/36	
%                          0       1/72	 1/36	 0       1/36	 1/18],...
%                          36,1);
%  
%            
%         C3 = reshape([  -1/24	-1/48	-1/48	 0       0       0	
%                         -1/48	-1/24	-1/48	 0       0       0	
%                         -1/48	-1/48	-1/24	 0       0       0	
%                          0       0       0       1/24	 1/48	 1/48	
%                          0       0       0       1/48	 1/24	 1/48	
%                          0       0       0       1/48	 1/48	 1/24 ],...
%                           36,1);

        C1 = reshape([      -8     8     0    -4     4     0
                            -8     8     0    -4     4     0
                            -8     8     0    -4     4     0
                            -4     4     0    -8     8     0
                            -4     4     0    -8     8     0
                            -4     4     0    -8     8     0] / 144,...
                        36,1);
                  
         
        C2  = reshape([     -8     0     8    -4     0     4
                            -8     0     8    -4     0     4
                            -8     0     8    -4     0     4
                            -4     0     4    -8     0     8
                            -4     0     4    -8     0     8
                            -4     0     4    -8     0     8]/144,...
                         36,1);
 
           
        C3 = reshape([      -6    -3    -3     6     3     3
                            -3    -6    -3     3     6     3
                            -3    -3    -6     3     3     6
                            -6    -3    -3     6     3     3
                            -3    -6    -3     3     6     3
                            -3    -3    -6     3     3     6 ]/144,...
                         36,1);


        
        % the element masss matrix
        M = reshape([ 1/36	 1/72	 1/72	 1/72	 1/144	 1/144	
                        1/72	 1/36	 1/72	 1/144	 1/72	 1/144	
                        1/72	 1/72	 1/36	 1/144	 1/144	 1/72	
                        1/72	 1/144	 1/144	 1/36	 1/72	 1/72	
                        1/144	 1/72	 1/144	 1/72	 1/36	 1/72	
                        1/144	 1/144	 1/72	 1/72	 1/72	 1/36  ],...
                        36,1);
        
        % the RHS element vector        
        F = 1/12 * [1
                    1
                    1
                    1
                    1
                    1];     
       
        
        % how to select points in gridObj.t
        idx = 1:6;
    end  
    
    methods(Access=public)
        function [idxvec0,idxvec1,idxvec2] = makeIndex(obj,idx,np)
            % helper function to making the sparse magic: re-aranging the
            % index of the matix...depends on the
            idxvec0 = reshape(idx,1,np*sqrt(size(obj.S1,1)));
            idxvec1 = reshape([idx;idx;idx;idx;idx;idx],1,np*size(obj.S1,1));            
            idxvec2 = reshape(...
                [idx(1,:);idx(1,:);idx(1,:);idx(1,:);idx(1,:);idx(1,:);...
                idx(2,:);idx(2,:);idx(2,:);idx(2,:);idx(2,:);idx(2,:);...
                idx(3,:);idx(3,:);idx(3,:);idx(3,:);idx(3,:);idx(3,:);...
                idx(4,:);idx(4,:);idx(4,:);idx(4,:);idx(4,:);idx(4,:);...
                idx(5,:);idx(5,:);idx(5,:);idx(5,:);idx(5,:);idx(5,:);...
                idx(6,:);idx(6,:);idx(6,:);idx(6,:);idx(6,:);idx(6,:)],...
                1,size(obj.S1,1)*np);   
        end  
        
        
        
        
    end
    
    methods(Static,Access=public)
        function [Q,G,H,R] = assemb(gridObj)
            % [Q,G,H,R] = assemb(gridObj)
            % Static method for bilinear 3D elemnts
            % assembles the boundary condition matrix
            %
            % (c) 2013 by Uwe Prüfert
            if isempty(gridObj.b)
                MException('bilinear3D:emptyboundary',...
                    'The boundary field is empty, initalize it!').throw;
            end

            nBoundaryElements = gridObj.nEdges;
            nTri = sum(gridObj.e(4,:)==0);
            nSqr = sum(gridObj.e(4,:)>0);
            nPts = size(gridObj.p,2);
            qval = zeros(1,nBoundaryElements);
            gval = zeros(1,nBoundaryElements);
            hval = zeros(1,nBoundaryElements);
            rval = zeros(1,nBoundaryElements);
 

            % Is this loop  nesseccary?
            for k = 1:nBoundaryElements     
                if gridObj.e(4,k)==0
                    [qval(k),gval(k),hval(k),rval(k)] = gridObj.bCoefficients(...
                        gridObj.p(:,gridObj.e(1:3,k)),gridObj.b(:,gridObj.e(5,k)));
                else
                    [qval(k),gval(k),hval(k),rval(k)] = gridObj.bCoefficients(...
                        gridObj.p(:,gridObj.e(1:4,k)),gridObj.b(:,gridObj.e(5,k)));
                end
            end

 
            % bdNum is the number of boundary segment...

            % all elements of boundary no bd
            % now check if we have triangles or squares
            % triangles:
            indxTri = find(gridObj.e(4,:)==0);

            % x-y-coordinates of the three face points... 
            % p = gridObj.p(1:2,gridObj.e(1:3,k));
            p1 = gridObj.p(1:2,gridObj.e(1,indxTri));
            p2 = gridObj.p(1:2,gridObj.e(2,indxTri));
            p3 = gridObj.p(1:2,gridObj.e(3,indxTri));

            % get element corner coordinates
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);    
            y31 = p3(2,:)-p1(2,:);

            J = x21.*y31-x31.*y21;


            Qe = reshape(lagrange12D.M*(J.*qval(indxTri)),1,9*nTri);
            ge = reshape(lagrange12D.F*(J.*gval(indxTri)),1,3*nTri); 

            He = reshape([1 0 0 0 1 0 0 0 1]'*(J.*hval(indxTri)),1,9*nTri);  
            re = reshape([1 1 1]'*(J.*rval(indxTri)),1,3*nTri); 

            % rearanging like in assema...
            indxTriPts = gridObj.e(1:3,indxTri);
            % 
            [indx0,indx1,indx2] = lagrange12D.makeIndex(indxTriPts,nTri);        

            Q = sparse(indx1,indx2,Qe,nPts,nPts); 
            G = sparse(indx0,1,ge,nPts,1);           
            H = sparse(indx1,indx2,He,nPts,nPts);
            R = sparse(indx0,1,re,nPts,1); 

            % Case 2. SQUARES
            indxSqr = find(gridObj.e(4,:)>0);

            indxSqrPts = gridObj.e(1:4,indxSqr);

            pSqr1 = gridObj.p(:,gridObj.e(1,:));
            pSqr2 = gridObj.p(:,gridObj.e(2,:));
            pSqr3 = gridObj.p(:,gridObj.e(3,:));

            % we set p1 = [0,0]
            p2 = [sqrt((pSqr2(1,indxSqr)-pSqr1(1,indxSqr)).^2+(pSqr2(2,indxSqr)-pSqr1(2,indxSqr)).^2)
                    zeros(size(indxSqr))];
            p3 = [zeros(size(indxSqr))
                 pSqr3(3,indxSqr)-pSqr2(3,indxSqr)];

             h = [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]';   
             r = [1 1 1 1]'; 

            J = p2(1,:).*p3(2,:)-p3(1,:).*p2(2,:);
      
            Qe = reshape((bilinear2D.M1+bilinear2D.M2)*(J.*qval(indxSqr)),1,16*nSqr);
            ge = reshape((bilinear2D.F1+bilinear2D.F2)*(J.*gval(indxSqr)),1,4*nSqr);  
            He = reshape(h*(J.*hval(indxSqr)),1,16*nSqr); 
            re = reshape(r*(J.*rval(indxSqr)),1,4*nSqr); 
 
            [indx0,indx1,indx2] = makeIndex(indxSqrPts,nSqr);

            Q = Q + sparse(indx1,indx2,Qe,nPts,nPts); 
            G = G + sparse(indx0,1,ge,nPts,1);           
            H = H + sparse(indx1,indx2,He,nPts,nPts);
            R = R + sparse(indx0,1,re,nPts,1);   

 
            % local function
                function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
                    idxvec0 = reshape(idx,1,np*4);
                    idxvec1 = reshape([idx;idx;idx;idx],1,np*16);
                    idxvec2 =  reshape([idx(1,:);idx(1,:);idx(1,:);idx(1,:);...
                        idx(2,:); idx(2,:); idx(2,:);idx(2,:);...
                        idx(3,:);idx(3,:);idx(3,:);idx(3,:);...
                        idx(4,:);idx(4,:);idx(4,:);idx(4,:)],1,16*np);     
                end 

        end 
    end
end
