function varargout = crossentropy_internal(varargin)
% CROSSENTROPY
%
% y = CROSSENTROPY(x,y)
%
% Computes/declares cross entropy -sum(x.*log(y))
%
% See also ENTROPY, KULLBACKLEIBLER

switch class(varargin{1})

    case 'double'    
        if nargin == 1
            % YALMIP flattens internally to [x(:);y(:)]
            z = varargin{1};
            z = reshape(z,[],2);
            x = z(:,1);
            y = z(:,2);
        else
            x = varargin{1}(:);
            y = varargin{2}(:);
        end
        ce = x(:).*log(y(:));
        ce(x==0) = 0;
        ce = real(ce);
        varargout{1} = -sum(ce);
        
    case 'char'
               
        operator = CreateBasicOperator('callback');
        operator.range = [-inf inf];
        operator.domain = [0 inf];       
        operator.derivative = @derivative;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function df = derivative(x)

z = reshape(x,[],2);
x = z(:,1);
y = z(:,2);

df = [-log(y);-x./y];


