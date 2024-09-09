classdef gridMethods3D
    % Class that implements some methods for  tetrahedrals
    % It should be used as an additional parent for class grid3D
    % Implements refine
    % refine 
    % 
    
    properties
        elementgeneration % marke refined elements
    end
    
    methods(Access = private)
        function [obj]=refine(obj,markedtri)
            % typical call:
            % [elem,egen,coord,diri,neum]=...
            %    refine(elem,egen,coord,diri,neum,marked)
            % marked=elem-numbers.
            % kemet  P-E-T
            % tri = t
            % po = p
            % varargin{1} = e
            % varargin{end} = market

                nB=1;
                p = obj.p(1:4,:)';
                e = obj.e'; % Noch varargin{?}
                t = obj.t';
                markedtri=unique(markedtri);
                
                %*** Wrapper−Schleife
                while ~isempty(markedtri)
                    marked=markedtri(1);
                    nE=size(t,1);
                    element2neighbors=provideNeighbors(t);
                    refined=[];
                    K=[];
                    F  =marked;
                    while (~isempty(F))
                        Fnew=[];
                        for T1=F
                            for T2=setdiff(element2neighbors(T1,[2,3]),[F,K,0])%** T2 ∈ N (P, T1 )\(F ∪K)
                                if (gen(T2)==gen(T1))
                                    Fnew=union(Fnew,T2);%** kompatibel teilbar
                                else
                                    [t,gen,p,varargout{1:nB},children,refinedtri]=...
                                        refine(t,gen,p,varargout{1:nB},T2);
                                    refined=[refined,refinedtri];
                                    element2neighbors=provideNeighbors(t);
                                    nE=size(t,1);
                                    if length(intersect(t(T1,:),t(children(1),:)))==3
                                        Fnew=union(Fnew,children(1));%** linkes Kind
                                    else
                                        Fnew=union(Fnew,children(2));%** rechtes Kind
                                    end
                                end
                            end
                        end
                        K=union(K,F); F=Fnew;
                    end
                    %*** Erzeugung des neuen Knotens
                    p(end+1,:)=1/2*( p(t(K(1),1),:) + p(t(K(1),4),:)); 
                    newNode=size(p,1);
                    %*** Verfeinerung der zur Kante gehoerigen Elemente
                    K=[marked,setdiff(K,marked)];%*** Damit children==[marked,nE+1] gilt
                    if (mod(gen(K),3)==0)
                        t([K,(nE+(1:length(K)))],:)=[ t(K,1),newNode*ones(length(K),1),t(K,[2,3]);
                        t(K,4),newNode*ones(length(K),1),t(K,[3,2]) ];
                    else
                        t([K,(nE+(1:length(K)))],:) = ...
                             [ t(K,1),newNode*ones(length(K),1),t(K,[2,3]);...
                        t(K,4),newNode*ones(length(K),1),t(K,[2,3]) ];
                    end
                    gen([K,nE+(1:length(K))])=[ gen(K)+1,gen(K)+1 ];
            %*** Verfeinerung der zur Kante gehoerigen Flaechen
                    if 0
                        for j=1:min(nargout-3,nB)
                            if ~isempty(varargout{j})
                                nBE=size(varargout{j},1); 
                                markedFaceNode=( (varargout{j}==tri(K(1),1)) | (varargout{j}==tri(K(1),4)) );
                                bisec=find(markedFaceNode(:,1) & markedFaceNode(:,3));
                                if ~isempty(bisec)
                                    varargout{j}([bisec,nBE+(1:length(bisec))'],:)=...
                                    [ varargout{j}(bisec,1),newNode*ones(length(bisec),1),...
                                      varargout{j}(bisec,2);...
                                      varargout{j}(bisec,3),newNode*ones(length(bisec),1),...
                                      varargout{j}(bisec,2) ];
                                end
                            end
                        end
                    end
                    %*** Hintergrundinformationen fuer die Rekursion
                    children=[marked,nE+1]; 
                    varargout{nB+1}=children; 
                    refined=[refined,K];
                    varargout{nB+2}=refined; 
                    markedtri=setdiff(markedtri,refined);
                end
                tr=TriRep(t,p); 
                varargout{1}=tr.freeBoundary; 
            end

        
    end
    methods(Access = private)
        
    end
end

