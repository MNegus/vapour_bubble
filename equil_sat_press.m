function sat_press = equil_sat_press(rad, params)
%EQUIL_SAT_PRESS Saturation pressure at equilibrium
%   By equilibrium, we mean no pertubations of temperature, but we can
%   still vary radius.

sat_press = params.refpress * exp((params.latent / params.Rspec) ...
    * (1 / params.reftemp - 1 / params.inftemp) ...
    - 2 * params.surftens / (params.Rspec * params.lden * params.inftemp * params.rad_0 * rad));

end

