function params_struct = plesset_zwick_params()
%PLESSET_ZWICK_PARAMS Physical parameters to compare to Plesset & Zwick
%   Constants independent of numerical scheme. Chosen the 10 degrees case.
%   Values taken on the saturation curve for the given temperature.

params_struct = struct;

% Reference temperature and pressure for use in Clausius
params_struct.reftemp = 373.15;
params_struct.refpress = 101325;

params_struct.Rspec = 461.52; % Specific gas constant

params_struct.latent = 2257e3; % Latent heat of evaporation

params_struct.inftemp = 375.15; % Temperature at infinity
params_struct.infpress = 101325; % Atmospheric pressure
params_struct.surftens = 0.058525; % Surface tension

params_struct.evap = 1; % Evaporation coefficient
params_struct.conden = 1; % Condensation coefficient

% Liquid parameters
params_struct.lden = 956.90; % Liquid density
params_struct.lvisc = 0.00027593; % Liquid viscosity
params_struct.ltherm = 0.679711; % Liquid thermal conductivity
params_struct.lspec = 4218.0; % Liquid specific heat capacity

% Vapour parameters
params_struct.vvisc = 1.2338e-05; % Vapour viscosity
params_struct.vbulkvisc = 6.5 * params_struct.vvisc; % Bulk viscosity (M.S. Cramer)
params_struct.vtherm = 0.025320; % Vapour thermal conductivity
params_struct.vspec = 2088.3; % Vapour specific heat capacity

% Calculating equilibrium radius
    function val =  equil_eqn(rad_0)
        % Equation for equilibrium radius
        val = 2 * params_struct.surftens / rad_0 + params_struct.infpress ...
            - params_struct.refpress * exp((params_struct.latent / params_struct.Rspec) ...
            * (1 / params_struct.reftemp - 1 / params_struct.inftemp) ...
            - 2 * params_struct.surftens / (params_struct.Rspec * params_struct.lden * params_struct.inftemp * rad_0));
    end

options = optimoptions(@fsolve); % Solver options

fun = @(a) equil_eqn(a); % Used in fsolve
guess = 10^-3; % Approximate guess for radius
params_struct.rad_0 = fsolve(fun, guess, options); % Solve for radius

% More vapour parameters
params_struct.vrfpress = params_struct.refpress * exp((params_struct.latent / params_struct.Rspec) ...
            * (1 / params_struct.reftemp - 1 / params_struct.inftemp) ...
            - 2 * params_struct.surftens / (params_struct.Rspec * params_struct.lden * params_struct.inftemp));
params_struct.vrfden = params_struct.vrfpress / (params_struct.Rspec * params_struct.inftemp);

% Typical value approximations
params_struct.U = 1; % Typical velocity, let's assume 1 to start
params_struct.t_0 = 10^-3; % Typical time scale, assume a millisecond
end

