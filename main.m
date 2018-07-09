function main()
% Create structs for physical and numerical parameters and load them into
% workspace
phys_params = plesset_zwick_params();

params_struct = numerical_params(phys_params);

% Unpack required parameters
reftemp = params_struct.reftemp; % Reference temperature used in Clausius
refpress = params_struct.refpress; % Reference pressure used in Clausius

Rspec = params_struct.Rspec; % Specific gas constant
latent = params_struct.latent; % Latent heat of evaporation

inftemp = params_struct.inftemp; % Temperature at infinity
infpress = params_struct.infpress; % Atmospheric pressure
surftens = params_struct.surftens; % Surface tension

evap = params_struct.evap; % Evaporation coefficient
conden = params_struct.conden; % Condensation coefficient

% Liquid parameters
lden = params_struct.lden; % Liquid density
lvisc = params_struct.lvisc; % Liquid viscosity
ltherm = params_struct.ltherm; % Liquid thermal conductivity
lspec = params_struct.lspec; % Liquid specific heat capacity

% Vapour parameters
vvisc = params_struct.vvisc; % Vapour viscosity
vbulkvisc = params_struct.vbulkvisc; % Bulk viscosity (M.S. Cramer)
vtherm = params_struct.vtherm; % Vapour thermal conductivity
vspec = params_struct.vspec; % Vapour specific heat capacity
vrfpress = params_struct.vrfpress; % Vapour reference pressure
vrfden = params_struct.vrfden; % Vapour reference density

% Typical velocity and time
U = params_struct.U;
t_0 = params_struct.t_0;
rad_0 = params_struct.rad_0;

% Numerical parameters
eps = params_struct.eps; % The small epsilon parameter in equations
M = params_struct.M; % Number of nodes in vapour
N = params_struct.N; % Total number of spacial nodes
dx = params_struct.dx; % Grid size 
total_time = params_struct.total_time; % Total time to run

% Defines ranges in collected array for state variables
vden_range = 1 : M;
vvel_range = M + 1 : 2 * M;
vtemp_range = 2 * M + 1 : 3 * M;
ltemp_range = 3 * M + 1 : N + 2 * M + 1;
lboundvel_idx = N + 2 * M + 2;
rad_idx = N + 2 * M + 3;


% Creates collected arrays for initial conditions
X_0 = zeros(rad_idx, 1); % Collected array for variables
X_0(rad_idx) = 1; % Specify that the initial (dimensionless) radius is 1

Xp_0 = zeros(rad_idx, 1); % Collected array for time derivatives
Xp_0(rad_idx) = 0; % Initial growth speed of bubble


% Useful arrays
vap_bulk = rscale(2 : M - 1); % Spacial distances in bulk vapour
liq_bulk = rscale(M + 1: N - 1); % Spacial distances in bulk liquid



% Useful constants
vel_balance = U * t_0 / rad_0; % Velocity term appearing a lot

vvisc_term = t_0 * (4 * vvisc / 3 + vbulkvisc) ...
    / (vrfden * rad_0^2); % Viscosity term for vapour

vtherm_term = vtherm * t_0 / (rad_0^2 * vrfden * vspec); 
    % Thermal conduction term for vapour

ltherm_term = ltherm * t_0 / (rad_0^2 * lden * lspec);
    % Thermal conduction term for liquid

% Attempts to find proper initial condition
fixed_X_0 = zeros(rad_idx, 1);
fixed_X_0(rad_idx) = 1;

fixed_Xp_0 = zeros(rad_idx, 1);
fixed_Xp_0(rad_idx) = 1;

% Compute accurate initial conditions
% [X_0, Xp_0] = decic(@odefun, 0, X_0, fixed_X_0, Xp_0, fixed_Xp_0);

save('X_0.mat', 'X_0');
save('Xp_0.mat', 'Xp_0');

options = odeset('Stats','on', 'RelTol', 1e-4);
tspan = [0, 1];
[t, y] = ode15i(@odefun, tspan, X_0, Xp_0, options);
solution = [t, y];

save('vars.mat', 'solution');
save('time.mat', 't');

plot(t, y(:, rad_idx));

    function res = odefun(t, X, Xp)
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining arrays
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        res = ones(rad_idx, 1); % Residues array
        
        % Creating arrays for state variables
        vden = X(vden_range).'; % Vapour density
        vden_deriv = Xp(vden_range).';
        
        vvel = X(vvel_range).'; % Vapour velocity
        vvel_deriv = Xp(vvel_range).';
        
        vtemp = X(vtemp_range).'; % Vapour temperature
        vtemp_deriv = Xp(vtemp_range).';
        
        ltemp = X(ltemp_range).'; % Liquid temperature
        ltemp_deriv = Xp(ltemp_range).';
        
        lboundvel = X(lboundvel_idx); % Liquid velocity at boundary
        lboundvel_deriv = Xp(lboundvel_idx);
        
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
        % Conservation of mass
        res(vden_range(2 : M - 1)) = vap_bulk .* vden_deriv(2 : M - 1) ...
            - (rad_deriv / rad) * (vap_bulk.^2) .* central_diff(vden, dx) ...
            + (vel_balance / rad) * (vap_bulk .* central_diff(vvel, dx) + 2 * vvel(2 : M - 1));
        
        % Conservation of momentum
        res(vvel_range(2 : M - 1)) = vap_bulk.^2 .* vvel_deriv(2 : M - 1) ...
            - (vvisc_term / rad^2) * vap_bulk.^2 .* second_diff(vvel, dx) ...
            - ((rad_deriv / rad) * vap_bulk.^3 + (2 * vvisc_term / rad^2) .* vap_bulk) ...
                .* central_diff(vvel, dx) ...
            - (2 * vvisc_term / rad^2) * vvel(2 : M - 1) ...
            + t_0 * vrfpress / (rad_0 * U * vrfden) ...
                * (1 / rad) * vap_bulk.^2 .* central_diff(vpress(vden, vtemp), dx);
        
        % Conservation of energy
        res(vtemp_range(2 : M - 1)) = vap_bulk .* vtemp_deriv(2 : M - 1) ...
            - ((rad_deriv / rad) * vap_bulk.^2 + (2 * vtherm_term / rad^2)) ...
                .* central_diff(vtemp, dx) ...
            - (vtherm_term / rad^2) * vap_bulk .* second_diff(vtemp, dx) ...
            + vrfpress * U * t_0  / (rad_0 * vrfden * vspec * inftemp) ...
                * (1 / rad) * (vap_bulk .* central_diff(vvel, dx) + 2 * vvel(2 : M - 1));
     
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Boundary conditions at bubble
        % %%%%%%%%%%%%%%%%%%%%%%%%%%   
        
        % First conservation of mass
        res(vden_range(M)) = t_0 / (vrfden * rad_0) ...
            * massflux(rad, ltemp(1), vtemp(M), vpress(vden(M), vtemp(M))) ...
            + eps * vel_balance * vvel(M) - (1 + eps * vden(M)) * rad_deriv;
            
        % Second conservation of mass
        res(rad_idx) = t_0 / (lden * rad_0) ...
            * massflux(rad, ltemp(1), vtemp(M), vpress(vden(M), vtemp(M))) ...
            + eps * vel_balance * lboundvel - rad_deriv;
        
        % Conservation of momentum, used to calculate liquid pressure at boundary
        lboundpress = vrfpress / infpress ...
            - U * lden * rad_0 / (t_0 * infpress) * (lboundvel - vvel(M)) * rad_deriv ...
            - U / (infpress * rad_0 * rad) * (4 * lvisc * lboundvel ...
                    + (4 * vvisc / 3) * ((vvel(M) - vvel(M-1)) / dx - vvel(M)) ...
                    + vbulkvisc * ((vvel(M) - vvel(M-1)) / dx + 2 * vvel(M))) ...
            - (1 / eps) * (1 - vrfpress / infpress  ...
                + 2 * surftens / (infpress * rad_0 * rad));
        
        
%         res(lboundpress_idx) = 1 - vrfpress / infpress ...
%             + 2 * surftens / (infpress * rad_0 * rad) ...
%             + eps * (lboundpress - vrfpress / infpress ...
%                 + U * lden * rad_0 / (t_0 * infpress) * (lboundvel - vvel(M)) * rad_deriv ...
%                 + U / (infpress * rad_0 * rad) * (4 * lvisc * lboundvel ...
%                     + (4 * vvisc / 3) * ((vvel(M) - vvel(M-1)) / dx - vvel(M)) ...
%                     + vbulkvisc * ((vvel(M) - vvel(M-1)) / dx + 2 * vvel(M))));
        
        % Conservation of energy
        res(vtemp_range(M)) = (2 * surftens / (rad_0 * lden * lspec * inftemp * rad) ...
            - (1 - vspec / lspec)) * rad_deriv ...
            + eps * (vel_balance * (1 - vspec / lspec) * lboundvel ...
                - (ltemp(1) - vspec * vtemp(M) / lspec) * rad_deriv ...
                + U * t_0 / (rad_0 * lden * lspec * inftemp) * (infpress * lboundvel ...
                    - vrfpress * vvel(M)) ...
                - t_0 / (rad_0^2 * lden * lspec * rad) ...
                    * (ltherm * (ltemp(2) - ltemp(1)) / dx ...
                        - vtherm * (vtemp(M) - vtemp(M-1)) / dx));
        
        % Temperature balance
        res(ltemp_range(1)) = ltemp(1) - vtemp(M);
        
        % Liquid boundary velocity
        res(lboundvel_idx) = lboundpress - (lden / infpress) ...
            * ((U * rad_0 / t_0) * (lboundvel_deriv * rad + 2 * lboundvel * rad_deriv) ...
                - eps * U^2 * lboundvel^2 / 2);
        
        % Density bulk equation applied at the boundary
        res(vvel_range(M)) = vden_deriv(M) ...
            - (rad_deriv / rad) * (vden(M) - vden(M - 1)) / dx ...
            + vel_balance / rad * ((vvel(M) - vvel(M - 1)) / dx + 2 * vvel(M));
        
%         % Liquid temperature bulk equation at the boundary
%         res(vvel_range(M)) = ltemp_deriv(1) ...
%             - (rad_deriv / rad + 2 * ltherm_term / rad^2) * (ltemp(2) - ltemp(1)) / dx ...
%             - ltherm_term / rad^2 * (ltemp(3) - 2 * ltemp(2) + ltemp(1)) / dx^2;
                    
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Liquid bulk
        % %%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        % Conservation of energy
        res(ltemp_range(2 : N - M)) = ltemp_deriv(2 : N - M) ...
            - ((rad_deriv / rad) * liq_bulk + (2 * ltherm_term / rad^2) ./ liq_bulk) ...
                .* central_diff(ltemp, dx) ...
            - (ltherm_term / rad^2) * second_diff(ltemp, dx);
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Boundary at infinity
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Temperature pertubation is zero at infinity
        res(ltemp_range(N - M + 1)) = ltemp(N - M + 1);
        
        
        
    end

    % Ideal gas law for pressure in vapour
    function val = vpress(vden, vtemp)
        val = vden + vtemp;
    end

    % Saturation pressure at equilibrium
    function sat_press = equil_sat_press(rad)
    sat_press = refpress * exp((latent / Rspec) ...
    * (1 / reftemp - 1 / inftemp) ...
    - 2 * surftens / (Rspec * lden * inftemp * rad_0 * rad));
    end


    % Mass flux from HKS equation
    function j = massflux(rad, lboundtemp, vboundtemp, vboundpress)
        j = 2 / (2 - conden) * sqrt(1 / (2 * pi * inftemp * Rspec)) ...
            * (evap * equil_sat_press(rad) - conden * vrfpress ...
            + eps * (evap * equil_sat_press(rad) * lboundtemp / (Rspec * inftemp) ...
                * (latent - 2 * surftens / (lden * rad_0 * rad) - Rspec * inftemp / 2) ...
                - conden * vrfpress * (vboundpress - vboundtemp / 2)));
    end

    % Creates function for the j^th rscale (physical position) 
    function position = rscale(j)
        position = dx * (j-1);
    end

    

end