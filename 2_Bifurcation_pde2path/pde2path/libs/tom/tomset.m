function options = tomset(varargin)
% tomset: from TOM
% options = tomset(varargin)
%TOMSET  Create/alter TOM OPTIONS structure.
%   OPTIONS = TOMSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is 
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names. 
%   
%   OPTIONS = TOMSET(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS. 
%   
%   OPTIONS = TOMSET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties overwrite 
%   corresponding old properties. 
%   
%   TOMSET with no input arguments displays all property names and their
%   possible values. 
%   
%TOMSET PROPERTIES
%   
%RelTol - Relative tolerance for the error/residual [  scalar {1e-3} ]
%   On eachpoint  of the mesh,
%   component i of the error satisfies  
%          norm( err(i) / max( RelTol* [abs(y(i)) , AbsTol(i)/RelTol] ) ) <= 1.
%
%Abstol - Absolute tolerance for the error/residual [ vector or scalar {1e-6} ]
%   This scalar applies to all components of the error/residual vector, and
%   defaults to 1e-6.   On eachpoint  of the mesh,
%   component i of the error satisfies  
%          norm( err(i) / max( RelTol* [abs(y(i)) , AbsTol(i)/RelTol] ) ) <= 1.
%
%Nmax - Maximum number of mesh points allowed [positive integer {floor(1000/n)}]
%
%Stats - Display final computational cost statistics  [ on | {off} ]
%
%Stats_step - Display  information at each step statistics  [ on | {off} ]
%
%PrintG - Plot approximate solution at each step  [ on | {off} ]
%
%IndexG - Index of the solution to be plotted [ {1} | 2 | ... | n ]
%
%Itnlmax - maximum number of non linear iteration [positive integer {50}]
%
%Itlinmax - maximum number of linear iteration [positive integer {50}]
%
%Vectorized - Vectorized ODE function  [ on | {off} ]
%   Set this property 'on' if the derivative function 
%   ODEFUN([x1 x2 ...],[y1 y2 ...]) returns [ODEFUN(x1,y1) ODEFUN(x2,y2) ...].  
%   When parameters are present, the derivative function
%   ODEFUN([x1 x2 ...],[y1 y2 ...],p) should return 
%   [ODEFUN(x1,y1,p) ODEFUN(x2,y2,p) ...].  
%
%Monitor - mesh selection [ {1} | 2| 3]
% monitor == 1 conditioning and error
% monitor == 2 approx conditioning and error
% monitor == 3 error
%
%Stabcondpar - on force the conditioning parameter to be stabilezed, used
% with monitor = 1 or monitor = 1
%
%
%   See also TOMGET,  TOM.
%
% Francesca Mazzia
% Modification of the bvpset file by:
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2003/06/23 17:12:47 $


% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('          RelTol: [ positive scalar  {1e-3} ]\n');
  fprintf('          AbsTol: [ positive scalar or vector {1e-6} ]\n');
  fprintf('       FJacobian: [ function ]\n');
  fprintf('      BCJacobian: [ function ]\n');
  fprintf('        ForceJAC: [ {on} | off ]\n');
  fprintf('           Order: [ 2 | {6}  ]\n');
  fprintf('           Stats: [ on | {off} ]\n'); 
  fprintf('      Stats_step: [ on | {off} ]\n');
  fprintf('          PrintG: [ on | {off} ]\n');
  fprintf('          IndexG: [ {1} | 2 | ... | n ]\n');
  fprintf('            Nmax: [ nonnegative integer {floor(1000/n)} ]\n'); 
  fprintf('         Itnlmax: [ nonnegative integer {50}]\n');
  fprintf('        Itlinmax: [ nonnegative integer {50}]\n');
  fprintf('      Vectorized: [ on | {off} ]\n'); 
  fprintf('         Monitor: [ {1} | 2 | 3]\n'); 
  fprintf('     Stabcondpar: [ on | {off} ]\n');
  fprintf('\n');
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
    'Itnlmax    '
    'Itlinmax   '
    'Nmax       '
    'Vectorized '
    'Monitor    '
    'Stabcondpar'
    ];
    
[m,n] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in tomset(o1,o2,...).
options = [];
i = 1;
while i <= nargin
  arg = varargin{i};
  if isstr(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(sprintf(['Expected argument %d to be a string property name '...
                     'or an options structure\ncreated with TOMSET.'], i));
    end
    if isempty(options)
      options = arg;
    else
      for j = 1:m
        val = getfield(arg,Names(j,:));
        if ~isequal(val,[])             % empty strings '' do overwrite
          options = setfield(options,Names(j,:),val);
        end
      end
    end
  end
  i = i + 1;
end
if isempty(options)
  for j = 1:m
    options = setfield(options,Names(j,:),[]);
  end
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;                      % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};    
  if ~expectval
    if ~isstr(arg)
      error(...
        sprintf('Expected argument %d to be a string property name.', i));
    end
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(sprintf('Unrecognized property name ''%s''.', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;                      % we expect a value next    
  else
    options = setfield(options,Names(j,:),arg);
    expectval = 0;      
  end
  i = i + 1;
end

if expectval
  error(sprintf('Expected value for property ''%s''.', arg));
end
