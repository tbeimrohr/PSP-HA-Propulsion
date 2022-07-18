function [mass] = finMass(density,rootChord,tipChord, span, thickness)
% finMass Calculates mass of a rocket fin
% This function caculates the mass of a trapesodial shaped fin.
%% -=-Inputs-=-
% 1. denisty: density of the material of the fin.
% 2. rootChord: length of the edge of the fin connected to the airframe.
% 3. tipChord: length of the outer edge of the fin.
% 4. span: distance from root chord to edge.
% 5. thickness: thickness of the fin

%% Calculations
area = ((rootChord + tipChord) / 2) * span;
volume = area * thickness;
mass = volume * density;
end
