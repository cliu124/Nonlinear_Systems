function [p,e,t]=hudistmesh(fd, fh, h0, box, iteration_max, pfix, varargin)
% local copy by HU, with small mods, of distmesh by Per-Olof Persson
    %  Parameters:
    %
    %    Input, function_handle FD  signed distance function d(x,y).
    %           function_handle FH  scaled edge length function h(x,y).
    %           double              H0, the initial smallest edge length.
    %           double              BOX(2,3), a bounding box for the region.
    %           double              ITERATION_MAX, the maximum number of iterations.
    %                               The iteration might terminate sooner than this
    %                               limit, if the program decides that the mesh 
    %                               has converged.
    %           double PFIX(NFIX,3) the coordinates of nodes that are
    %                               required to be included in the mesh.
    %           VARARGIN,           additional parameters that can be passed to FD
    %
    %    Output, double P(N,3)      node coordinates.
    %            double T(NT,4)     tetrahedron indices.
    ptol=0.001; ttol=0.1;  L0mult=1.1;  deltat=0.1; 

    geps=0.01*h0; deps=sqrt(eps)*h0; iteration=0;
    iteration_max=max(iteration_max,1);        
    cbox=cell(1,3); 
    for k=1:3;  cbox{k}=box(1,k):h0:box(2,k);  end
    pp=cell(1,3);
    [pp{:}]=ndgrid(cbox{:});
    p=zeros(numel(pp{1}),3);
    for k=1:3; p(:,k)=pp{k}(:); end
    % 2. Remove points outside the region, apply the rejection method.
    %
    p=p(feval(fd,p,varargin{:}) < geps,:);  r0=feval(fh,p);
    p=p(rand(size(p,1),1)<min(r0)^3./r0.^3,:); 
    p=[pfix;p]; N=size(p,1); 

    p0=inf;
    while(iteration < iteration_max) 
        iteration=iteration+1;
        if(ttol*h0<max(sqrt(sum((p-p0).^2,2))))
            p0=p; t=delaunayn(p); 
            pmid=zeros(size(t,1),3);
            for k=1:3+1; pmid=pmid+p(t(:,k),:) / (3+1); end
            %
            %  Only keep those simplices in the new triangulation for which
            %  the centroids is inside the region, or not too far outside.
            %
            t=t(feval(fd, pmid, varargin{:})< -geps,:);
            %
            %  4. Describe each edge by a unique pair of nodes.
            %
            pair=zeros(0,2); localpairs=nchoosek(1:3+1,2);
            for k=1:size(localpairs,1)
                pair=[pair;t(:,localpairs(k,:))];
            end
            pair=unique(sort(pair,2),'rows'); 
        end

        bars=p(pair(:,1),:)-p(pair(:,2),:);
        L=sqrt(sum(bars.^2,2));
        L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
        L0=L0*L0mult*(sum(L.^3)/sum(L0.^3))^(1/3);
        F=max(L0-L,0);  Fbar=[bars,-bars].*repmat(F./L,1,2*3);
        dp=full(sparse(pair(:,[ones(1,3),2*ones(1,3)]), ...
                ones(size(pair,1),1)*[1:3,1:3], ...
                Fbar,N,3));
        dp(1:size(pfix,1),:)=0; 
        p=p+deltat*dp;

        for steps=1:2
            d=feval(fd, p, varargin{:}); 
            indx=( 0 < d);
            gradd=zeros(sum(indx), 3);

            for k=1:3
                a=zeros(1,3);
                a(k)=deps;
                d1x=feval(fd, p(indx,:)+ones(sum(indx),1)*a, varargin{:});
                gradd(:,k)=( d1x - d(indx))/ deps;
            end

            p(indx,:)=p(indx,:) - d(indx)*ones(1,3) .* gradd;

        end

        maxdp=max(deltat*sqrt(sum(dp(d < -geps,:).^2, 2)));

        if(maxdp < ptol*h0)
            break;
        end

    end
    size(p), size(pfix)
    p=[pfix; p];  
    t=delaunayn(p); % a posteriori add pfix again 
    e=surf2tri(p, t);
    p=p'; e=e'; t=t';

    function tri=surf2tri(p, t)               
        %
        faces=[  t(:,[1,2,3]);
                   t(:,[1,2,4]);
                   t(:,[1,3,4]);
                   t(:,[2,3,4])];
        node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
        faces=sort(faces,2);
        [~,iindx,jindx]=unique(faces,'rows');
        vec=histc(jindx,1:max(jindx));
        qx=find(vec==1);
        tri=faces(iindx(qx),:);
        node4=node4(iindx(qx));
        %
        % Orientation
        %
        v1=p(tri(:,2),:)-p(tri(:,1),:);
        v2=p(tri(:,3),:)-p(tri(:,1),:);
        v3=p(node4,:)-p(tri(:,1),:);
        iindx=find(dot(cross(v1,v2,2),v3,2)>0);
        tri(iindx,[2,3])=tri(iindx,[3,2]);
    end
end

function d=ddiff(d1,d2)
    % DDIFF returns the signed distance to a region that is the difference of two regions.
    d=max(d1, -d2);             
end

function d=dintersect(d1, d2)
    % DINTERSECT sets the signed distance to the intersection of two regions.
    d=max(d1,d2);         
end
