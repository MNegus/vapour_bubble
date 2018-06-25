function params_struct = plesset_zwick_params()
%WATER_PARAMS Physical parameters to compare to Plesset & Zwick
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

% Liquid parameters
params_struct.lden = 956.90; % Liquid density
params_struct.lvisc = 0.00027593; % Liquid viscosity
params_struct.ltherm = 0.679711; % Liquid thermal conductivity
params_struct.lspec = 4218.0; % Liquid specific heat capacity

% Vapour parameters
params_struct.vvisc = 1.2338e-05; % Vapour viscosity

params_struct.vtherm = 0.025320; % Vapour thermal conductivity
params_struct.vspec = 2088.3; % Vapour specific heat capacity

% Equilibrium radius
    function val =  equil_eqn(rad_0)
        % Equation for equilibrium radius
        val = 2 * params_struct.surftens / rad_0 + params_struct.infpress ...
            - params_struct.refpress * exp((params_struct.latent / params_struct.Rspec) ...
            * (1 / params_struct.reftemp - 1 / params_struct.inftemp) ...
            - 2 * params_struct.surftens / (params_struct.Rspec * params_struct.lden * params_struct.inftemp * rad_0));
    end

options = optimoptions(@fsolve, 'Display', 'iter', ...
'StepTolerance', 10^-20); % Solver options

fun = @(a) equil_eqn(a); % Used in fsolve

guess = 10^-3; % Approximate guess for radius

params_struct.rad_0 = fsolve(fun, guess, options);

end

