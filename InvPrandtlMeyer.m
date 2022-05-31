function [M] = InvPrandtlMeyer(gamma,nu)
% This function numerically inverts the PM function to find M

%% Inputs:
% gamma
% nu (rad)
%% Outputs:
% M

% Check to make sure nu is within the bounds of what is possible
nu_min = 0.0;
nu_max = pi*0.5*(sqrt((gamma+1.0)/(gamma-1.0))-1.0);
tolerance = 1.0E-8;
if nu < nu_min
    error("Error in InvPrandtlMeyer.  nu is less than 0.0")
end
if nu > nu_max
    error("Error in InvPrandtlMeyer.  nu is greater than max value")
end

if nu==0
    M=1.0;
    return
end
% Determine upper bound
M_upper = 2.0;
nu_upper = PrandtlMeyer(gamma,M_upper);
while nu > nu_upper
    M_upper = M_upper*1.5;
    nu_upper = PrandtlMeyer(gamma,M_upper);
end
% Set lower bound
M_lower = 1.0;
converged=false;

M_guess = 0.5*(M_upper+M_lower);
nu_guess = PrandtlMeyer(gamma,M_guess);
nitit = 1;
while abs((nu_guess-nu)/nu)> tolerance
    if nu_guess<nu
        M_lower = M_guess;
    else
        M_upper = M_guess;
    end
    M_guess = 0.5*(M_upper+M_lower);
    nu_guess = PrandtlMeyer(gamma,M_guess);
    nitit=nitit+1;
end
 M = M_guess;
 return
    
