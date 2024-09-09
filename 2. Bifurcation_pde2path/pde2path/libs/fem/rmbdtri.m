function in=rmbdtri(p,in)
% RMBDTRI: remove boundary triangles from refinement-list
%
%  in=rmbdtri(p,in)
% po=points, t=triangles, in=refinement-list
% Important settings: p.nc.bddistx, p.nc.bddisty
%
% See also stanparam
[po,t]=getpte(p); x=po(1,:); y=po(2,:); x1=min(x); x2=max(x); y1=min(y); y2=max(y);
dx=p.nc.bddistx; dy=p.nc.bddisty; % cut off to bdry
ni=size(in,2); keep=[]; bcp=p.sw.bcper;  
if size(in,1)>0; 
for j=1:ni % loop over triangle indizes
    for k=1:3 % loop over triangle-points
        z=po(:,t(k,in(j))); % k-th-point of triangle j
        if(0)
            if (bcp==2 || bcp==3)
                if (abs(z(2)-y1)<dy || abs(z(2)-y2)<dy) keep=[keep in(j)]; break; end
            end
            if (bcp==1 || bcp==3)
                if (abs(z(1)-x1)<dx || abs(z(1)-x2)<dx) keep=[keep in(j)]; break; end
            end
        end
        for m=bcp
            if(m==1)
                if(abs(z(1)-x1)<dx || abs(z(1)-x2)<dx) keep=[keep in(j)]; end
            elseif(m==2)
                if (abs(z(2)-y1)<dy || abs(z(2)-y2)<dy) keep=[keep in(j)]; end
            end
        end
    end
end
in=setdiff(in,keep);
end

