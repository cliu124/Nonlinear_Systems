function [FillOp,DropOp,nu]=getPerOp(p,varargin)
% GETPEROP: generate matrices for changing from Neumann to periodic bc.
% mod for Xcont if p.sw.orgper~=0; otherwise calls org-version (default)  
if(nargin>1) dir=varargin{1}; else dir=p.sw.bcper; end
if dir==0; FillOp=1;DropOp=1; nu=p.nu; return; end 
try; org=p.sw.orgper; catch; org=1; end; %org=1; 
if org; [FillOp,DropOp,nu]=getPerOporg(p,dir); return; end 
p.pdeo.grid.p=p.X'; poi=getpte(p); n=size(poi,1); 
if(length(dir)>n || max(dir)>n || min(dir) < 1)
    fprintf(' Error in direction for getPerOp !!! \n');
end
[FillOp,DropOp,nu]=getPerOp1dir(p,dir(1));
if(length(dir)>1 && length(dir)<=3)
    q=p; q.np=size(DropOp,1)/q.nc.neq;
    np_old=size(DropOp,2)/q.nc.neq;   
    try q.pdeo.grid.p=transpose(DropOp(1:q.np,1:np_old)*q.pdeo.grid.p');
    catch q.mesh.p=transpose(DropOp(1:q.np,1:np_old)*q.mesh.p'); end
     
    for k=dir(2:end)
        [tmp_fill, tmp_drop, tmp_nu]=getPerOp1dir(q,dir(k));        
        FillOp=FillOp*tmp_fill; DropOp=tmp_drop*DropOp; nu=tmp_nu;        
        q.np=size(tmp_drop,1)/p.nc.neq;
        np_old=size(tmp_drop,2)/p.nc.neq;        
      %drop the boundary points in q which have just been periodically identified     
       try q.pdeo.grid.p=transpose(tmp_drop(1:q.np,1:np_old)*q.pdeo.grid.p');
       catch q.mesh.p=transpose(tmp_drop(1:q.np,1:np_old)*q.mesh.p'); end    
    end
end
if isfield(p,'pdeo');  % restore points in q=p (obj. orient. weirdness!) 
   q.pdeo.grid.p=poi; 
end 

