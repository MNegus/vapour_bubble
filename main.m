function out = main()
% Create structs for physical and numerical parameters and load them into
% workspace
phys_params = plesset_zwick_params();

params = numerical_params(phys_params);

% Unpack required parameters for brevity
M = params.M;
N = params.N;
dx = params.dx;

% Defines ranges in collected array for state variables
vden_range = 1 : M;
vvel_range = M + 1 : 2 * M;
vtemp_range = 2 * M + 1 : 3 * M;
ltemp_range = 3 * M + 1 : N + 2 * M + 2;
lboundvel_idx = N + 2 * M + 3;
lboundpress_idx = N + 2 * M + 4;
rad_idx = N + 2 * M + 5;

% Creates collected arrays for initial conditions
X_0 = zeros(1, rad_idx); % Collected array for variables
X_0(rad_idx) = 1; % Specify that the initial (dimensionless) radius is 1

Xp_0 = zeros(1, rad_idx); % Collected array for time derivatives
Xp_0(rad_idx) = 10^-4; % Initial growth speed of bubble

% Useful arrays
vap_bulk = rscale(2 : M - 1); % Spacial distances in bulk vapour
liq_bulk = rscale(M + 1: N - 1); % Spacial distances in bulk liquid


% Useful constants
vvisc_comb = 4 * params.vvisc / 3 + params.vbulkvisc; % Combined viscosity
vel_balance = params.U * params.t_0 / params.rad_0;

thing = odefun(1, X_0, Xp_0);
    function res = odefun(t, X, Xp)
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining arrays
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        res = ones(1, rad_idx); % Residues array
        
        % Creating arrays for state variables
        vden = X(vden_range); % Vapour density
        vden_deriv = Xp(vden_range);
        
        vvel = X(vvel_range); % Vapour velocity
        vvel_deriv = Xp(vvel_range);
        
        vtemp = X(vtemp_range); % Vapour temperature
        vtemp_deriv = Xp(vtemp_range);
        
        ltemp = X(ltemp_range); % Liquid temperature
        ltemp_deriv = Xp(ltemp_range);
        
        lboundvel = X(lboundvel_idx); % Liquid velocity at boundary
        lboundvel_deriv = Xp(lboundvel_idx);
        
        lboundpress = X(lboundpress_idx); % Liquid pressure at boundary
        
        rad = X(rad_idx); % Radius of bubble
        rad_deriv = Xp(rad_idx);
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Boundary at origin
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Spacial derivatives are zero, for smoothness
        res(vden_range(1)) = vden(2) - vden(1);
        res(vvel_range(1)) = vvel(2) - vvel(1);
        res(vtemp_range(1)) = vtemp(2) - vtemp(1);
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Vapour bulk
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        res(vden_range(2 : M - 1)) = vden_deriv(2 : M - 1) ...
            - (rad_deriv / rad) * vap_bulk .* central_diff(vden, dx) ...
            + (vel_balance / rad) * (central_diff(vvel, dx) + 2 * (vvel(2 : M - 1)) ./ vap_bulk);
        
    
            
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Boundary at infinity
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Temperature pertubation is zero at infinity
        res(ltemp_range(N - M + 1)) = ltemp(N - M + 1);
    end

    % Creates function for the j^th rscale (physical position) 
    function position = rscale(j)
        position = dx * (j-1);
    end

end