function [FillOp,DropOp,nu]=getPerOpX(p,varargin)
% GETPEROPX: generate fill/drop matrices for manifolds X
% mod of getPerOp from standard pBCs 
%
%  [FillOp,DropOp,nu]=getPerOp(p)
%  [FillOp,DropOp,nu]=getPerOp(p,bcper,dir)

if(nargin>1) dir=varargin{1}; else dir=p.sw.bcper; end
if dir==0; FillOp=1;DropOp=1; nu=p.nu; return; end 
com=varargin{2}; %AM
poi=getpte(p); n=size(poi,1); 
if com==dir(1); [FillOp,DropOp, nu]=getPerOp1dirX(p,dir(1)); % wegen den "-" 
else [FillOp,DropOp,nu]=getPerOporg(p,dir(1));
end
if(length(dir)>1 && length(dir)<=3)
    q=p; q.np=size(DropOp,1)/q.nc.neq;
    np_old=size(DropOp,2)/q.nc.neq;
    %drop the boundary points in q which have just been periodically identified
    try q.pdeo.grid.p=transpose(DropOp(1:q.np,1:np_old)*q.pdeo.grid.p');
    catch q.mesh.p=transpose(DropOp(1:q.np,1:np_old)*q.mesh.p'); end
     
    for k=dir(2:end)
        if com==k; [tmp_fill,tmp_drop, tmp_nu]=getPerOp1dirX(q,dir(k)); 
        else 
        [tmp_fill,tmp_drop, tmp_nu]=getPerOporg(q,dir(k));
        end
        FillOp=FillOp*tmp_fill; DropOp=tmp_drop*DropOp;  nu=tmp_nu;        
        q.np=size(tmp_drop,1)/p.nc.neq;
        np_old=size(tmp_drop,2)/p.nc.neq;
        
      %drop the boundary points in q which have just been periodically identified
     % q.pdeo.grid, q.pdeo.grid.p=transpose(tmp_drop(1:q.np,1:np_old)*q.pdeo.grid.p'); 17, pause 
       try q.pdeo.grid.p=transpose(tmp_drop(1:q.np,1:np_old)*q.pdeo.grid.p');
       catch q.mesh.p=transpose(tmp_drop(1:q.np,1:np_old)*q.mesh.p'); end    
    end
end
if isfield(p,'pdeo');  % restore points in q=p (obj. orient. weirdness!) 
   q.pdeo.grid.p=poi; 
end 

