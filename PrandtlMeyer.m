function [nu] = PrandtlMeyer(gamma,M)
% This function determines the Prandtl-Meyer Angle of a given flow.

%% Inputs:
% gamma = Specific Heat Ratio of fluid
% M = Mach Number of flow
% -------------------------------------------------------------------------
%% Outputs:
% nu = Prandtl-Meyer Angle (rad)
% -------------------------------------------------------------------------
% Prandtl-Meyer Functio
nu=sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(M^2-1)))-...
    atan(sqrt(M^2-1));
end