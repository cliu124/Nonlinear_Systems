function o = tomget(options,name,default)
% tom: extract options for TOM, by Francesca Mazzia
%
% o=tomget(options,name,default)
%TOMGET  Get TOM OPTIONS parameters.
%   VAL = TOMGET(OPTIONS,'NAME') extracts the value of the named property
%   from integrator options structure OPTIONS, returning an empty matrix if
%   the property value is not specified in OPTIONS. It is sufficient to type
%   only the leading characters that uniquely identify the property. Case is
%   ignored for property names. [] is a valid OPTIONS argument. 
%   
%   VAL = TOMGET(OPTIONS,'NAME',DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified
%   in OPTIONS. For example 
%   
%       val = tomget(opts,'RelTol',1e-4);
%   
%   returns val = 1e-4 if the RelTol property is not specified in opts.
%   
%   See also TOMSET, TOM
%
% Francesca Mazzia
% Modification of the bvpget file by:
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2003/06/23 17:12:59 $

if nargin < 2
  error('Not enough input arguments.');
end
if nargin < 3
  default = [];
end

if ~isempty(options) & ~isa(options,'struct')
  error('First argument must be an options structure created with TOMSET.');
end

if isempty(options)
  o = default;
  return;
end

Names = [
    'RelTol     '
    'AbsTol     '
    'FJacobian  '
    'BCJacobian '
    'ForceJAC   '
    'Order      '
    'Stats      '
    'Stats_step '
    'PrintG     '
    'IndexG     '
    'Nmax       '
    'Itnlmax    '
    'Itlinmax   '
    'Vectorized '
    'Monitor    '
    'Stabcondpar'
    ];

[m,n] = size(Names);
names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names); 
if isempty(j)               % if no matches
  error(sprintf(['Unrecognized property name ''%s''.  ' ...
                 'See TOMSET for possibilities.'], name));
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', name);
    msg = [msg '(' deblank(Names(j(1),:))];
    for k = j(2:length(j))'
      msg = [msg ', ' deblank(Names(k,:))];
    end
    msg = sprintf('%s).', msg);
    error(msg);
  end
end
%options, 
if any(strcmp(fieldnames(options),deblank(Names(j,:))))
 %  o = getfield(options, Names(j,:));
  try; o = getfield(options, deblank(Names(j,:))); catch; o=[]; end % HU
  if isempty(o)
    o = default;
  end
else
  o = default;
end
%o=20 % HU 

