function [FillOp,DropOp,nu]=getPerOp(p,varargin)
% GETPEROP: generate matrices for changing from Neumann to periodic bc.
%  Geometry is assumed to be a cuboid, bc to be Neumann 
%
%  [FillOp,DropOp,nu]=getPerOp(p)
%  [FillOp,DropOp,nu]=getPerOp(p,bcper)
%
% "FillOp" for modifying the stiffness matrix (and the load vector) for 
% periodic BCs. "DropOp" removes redundant entries for
% identified boundaries, "nu" is reduced length of vector u of unknowns.
% 
% See also getPerOp1dir, box2per, filltrafo.

% * bcper is a vector of size at most n (spatial dimension of the problem)
% * bcper=1: the x-ends of the interval to be identified (periodicity in x)
% * bcper=2: the y-ends of the interval to be identified (periodicity in y)
% * bcper=3: the z-ends of the interval to be identified (periodicity in z)
% * bcper=[1 2]: periodicity in x  and y
% * bcper=[1 3]: periodicity in x  and z
% * bcper=[2 3]: periodicity in y  and z
% * bcper=[1 2 3]: periodicity in x, y, and z

%p.nu, p.np, pause 
if(nargin>1) dir=varargin{1}; else dir=p.sw.bcper; end
if dir==0; FillOp=1;DropOp=1; nu=p.nu; return; end 

poi=getpte(p); n=size(poi,1);  %n=spatial dimension

if(length(dir)>n || max(dir)>n || min(dir) < 1)
    fprintf(' Error in direction for getPerOp !!! \n');
end

[FillOp,DropOp,nu]=getPerOp1dirorg(p,dir(1));

if(length(dir)>1 && length(dir)<=3)
    q=p; q.np=size(DropOp,1)/q.nc.neq;
    np_old=size(DropOp,2)/q.nc.neq;
    %drop the boundary points in q which have just been periodically identified
    try q.pdeo.grid.p=transpose(DropOp(1:q.np,1:np_old)*q.pdeo.grid.p');
    catch q.mesh.p=transpose(DropOp(1:q.np,1:np_old)*q.mesh.p'); end
     
    for k=dir(2:end)
        [tmp_fill, tmp_drop, tmp_nu]=getPerOp1dirorg(q,dir(k));
        
        FillOp=FillOp*tmp_fill;
        DropOp=tmp_drop*DropOp;
        nu=tmp_nu;
        
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

