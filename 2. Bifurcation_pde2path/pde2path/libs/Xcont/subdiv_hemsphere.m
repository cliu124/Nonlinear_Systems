function [V,F] = subdiv_hemsphere(iters,varargin)
% subdiv_hemsphere: generate hemisphere by subdiv and projection
% 
% mod of subdivided_sphere  from gptool 
  % default values
  radius = 1;
  subdivision_method = 'upsample';
  base = '';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Base','Radius','SubdivisionMethod'},{'base','radius','subdivision_method'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  switch base
  case 'icosahedron'
    %[V,F] = icosahedron();
      V = [eye(3);-[1 0 0; 0 1 0]];
  F = [
     1     2     3
     1     3     5
     2     4     3
     3     4     5];
  case 'octahedron'
    [V,F] = octahedron();
  otherwise
    % Compute the 12 vertices
    phi = (1+sqrt(5))/2;  % Golden ratio
    V = [   0   1  phi 
            0  -1  phi 
            1 0 0 
            -1 0 0
            1  phi  0  
          -1  phi  0  
            1 -phi  0  
          -1 -phi  0  
            phi 0   1  
          -phi 0   1  
          0 1 0
          0 -1 0];
    % Scale to required radius
    V = V/(sqrt(1+phi^2));
    % Define the adjacency matrix
    F = [1  2  9
          1  9  5
          1  5  11
          1  6  10
          1  10 2
          2  7  9
          10 8  2
          2  12  7
          3 5 9
          3 7 9
          4 8 10
          4 6 10
          1  6  11
          2  8 12 ];
  end

  V = normalizerow(V);
  for iter = 1:iters
    switch subdivision_method
    case 'upsample'
      [V,F] = upsample(V,F);
    case 'loop'
      [V,F] = loop(V,F);
    case 'sqrt3'
      [V,F] = sqrt3(V,F);
    end
    V = normalizerow(V);

  end
  V = V*radius;

end

