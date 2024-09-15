function b = boundarymatrix(q,g,h,r,varargin)
% BOUNDARYMATRIX: returns a 1D boundary matrix 
%
%  b = boundarymatrix(q,g,h,r,varargin)
% b = boundarymatrix([arglist]) by Uwe Pruefert
% 
%  usuage:
% same BC on both boundary:
% b = boundarymatrix            homogenious Neumann BCs  
% b = boundarymatrix(g)         Neumann BCs
% b = boundarymatrix(q,g)       Robin BCs 
% b = boundarymatrix([],[],h,r) Dirichlet BC
%
%  different BC: e.g. Dirichlet on the left bound, Robin on the right
% boundary
% b = boundarymatrix([],[],h,r,q,g,[],[])
% q,g,h,r are STRINGS like '1' 'sin(s)' etc... but NO inlines
%(c) Uwe.Pruefert@tu-berlin.de
%
% See also polygong, gnbc
%
% legacy setup. In OOPDE, use, e.g. bc=grid.robinBC(1,1); grid.makeBoundaryMatrix(bc); 

switch nargin                  % handle the "minimalistic syntax" option
    case {0}
        q = '0';
        g = '0';
        h = [];
        r = [];
    case {1}  
        g =  q;
        q = '0';
        h = [];
        r = [];
    case {2}
        h = [];
        r = [];
    otherwise
        %
end    
optargs = length(varargin);
nobdrs = max(1,1+optargs/4);    % number of definition blocks, if only one,
                                % we assume the same BCs on both boundaries

b = cell(nobdrs);               % initialize b as cell array, later we
                                % will overwrite it
if ~mod(optargs,4)==0,
    error('argument list has wrong lenght)')
end
   
if isempty(h)&&isempty(r),
   b1 = 1;                      % only Neumann BCs
   b1(2,1) = 0;
   b1(3,1) = length(q);
   b1(4,1) = length(g);
   b1=[b1;double(q)';double(g)'];
elseif isempty(q)&&isempty(g)   % Dirichlet BCs, setting the Neumann part  
   b1 = 1;                      % to zero
   b1(2,1) = 1;
        b1(3,1) = 1;
        b1(4,1) = 1;
        b1(5,1) = length(h);
        b1(6,1) = length(r);
        b1(7,1) = 48;
        b1(8,1) = 48;
        b1 = [b1;double(h)';double(r)'];    
end
b{1} = b1;

% more than one block of BCs in the parameter list
l = 1;

for k = 2:nobdrs,
    q = varargin{l}; l = l+1;
    g = varargin{l}; l = l+1;
    h = varargin{l}; l = l+1;
    r = varargin{l}; l = l+1;
    if isempty(h)&&isempty(r),
        b1 = 1; % 
        b1(2,1) = 0;
        b1(3,1) = length(q);
        b1(4,1) = length(g);
        b1 = [b1;double(q)';double(g)']; 
        b{k} = b1;
    elseif isempty(q)&&isempty(g) % 
        b1 = 1;
        b1(2,1) = 1;
        b1(3,1) = 1;
        b1(4,1) = 1;
        b1(5,1) = length(h);
        b1(6,1) = length(r);
        b1(7,1) = 48;
        b1(8,1) = 48;
        b1 = [b1;double(h)';double(r)'];
        b{k} = b1;
    end
end

                                % contruct the boundary matrix, we have to
                                % care about the right dimension of ALL
                                % rows when assembling  the matrix
                                
if nobdrs==1,                   % the case of one BCs for all boundaries
    bb=[b{1},b{1}];
else
    maxlengthb = 0;
    for l = 1:k,                % compute the max lenght of all entries of b
        if length(b{l})>maxlengthb,
            maxlengthb = length(b{l});
        end
    end
    bb = zeros(maxlengthb,nobdrs);  % initalize the matrix
    for l = 1:k,
        bb(1:length(b{l}),l)=b{l};
    end
end
b=bb;                           % return b
end
