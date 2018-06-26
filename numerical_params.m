function num_params = numerical_params()
    num_params.eps = 10^-4; % The small epsilon parameter in equations
    
    % Spacial parameters
    num_params.M = 200; % Number of nodes in vapour
    num_params.N = num_params.M * 5; % Total number of nodes
    num_params.dx = 1 / (num_params.M - 1);
    
    % Temporal parameters
    num_params.total_time = 10;
end