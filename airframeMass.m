function mass = airframeMass(density, innerDiameter, thickness, length)
% airframeMass calculates the mass of an an airframe
% the mass of the ariframe is calculated using the material density, inner
% diamater, outer diamater and lenth of the tube. The variables definitions
% and relavant unts are listed below.
%% -=-Inputs-=-
% 1. density: material density 
% 2. innerDiamater: inner diamater of the airframe tube 
% 3. thickness: radial thickness of the airframe tube
% 4. length: length or height of the airframe tube

%% Calculations
    outer_diameter = 2 * thickness + innerDiameter;
    mass = pi * length / 4 * (power(outer_diameter, 2) - power(innerDiameter, 2)) * density;

end
