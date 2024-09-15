function bc = gnbc(pneq,varargin)
% GNBC: returns a boundary matrix for generalized Neumann BC 
%         edgewise of form      n*c tensor grad(u) + q*u = g,
% where q is an pneq by pneq matrix and g a vector with pneq components. 
%
%  bc= gnbc(pneq,varargin)
%
% * For even length varargin:
%   varargin is vector of edgenum pairs (q,g), with one pair per edge 
%   of the domain.
% * For odd length varargin and equal (q,g) on each edge:
%   varargin = (edgenum, q, g)
%
% Note: q = [[12;6789];[345;101112]] means q21=6789 in PDEtool GUI.
%
% See also pdetool, polygong, boundarymatrix
if(mod(length(varargin),2)==0) 
    edgenum = length(varargin)/2;
else
    edgenum = varargin{1};
    if(length(varargin)~=3) 
        error('argument list has wrong length')
    end
    q = varargin{2};
    g = varargin{3};
    for jj=1:edgenum
        varargin{2*jj-1} = q;
        varargin{2*jj}   = g;
    end
end

bb(1,:)= pneq*ones(edgenum,1); % # of components 
bb(2,:)= zeros(edgenum,1); % set # Dirichlet b.c. to zero

for jj=1:edgenum,
    q = varargin{2*jj-1};
    g = varargin{2*jj};
    % Add lengths of entries of q and g to bb

    qa=q(:); % vector of entries of matrix q
    for j=1:length(qa),
        qs=mat2str(double(qa(j)));
        bb(j+2,jj)=length(qs);
        qc = uint8(qs); % ascii codes of characters as vector
        bq{j} = qc; 
    end

    ga=g(:); 
    for j=1:length(ga),
        gs=mat2str(double(ga(j)));
        bb(j+2+length(qa),jj)=length(gs);
        gc = uint8(gs); % ascii codes of characters as vector
        bg{j} = gc(:);
    end

    % Add ascii codes of characters of entries to bb

    bcount = 0;
    for j=1:length(qa),
        bbq = bq{j};
        for k=1:length(bbq)
            bb(k+2+length(qa)+length(ga)+bcount,jj) = bbq(k);
        end
        bcount = bcount + length(bbq);
    end

    for j=1:length(ga),
        bbg = bg{j};
        for k=1:length(bbg)
            bb(k+2+length(qa)+length(ga)+bcount,jj) = bbg(k);
        end
        bcount = bcount + length(bbg);
    end

    bc=bb; % return bb
end
