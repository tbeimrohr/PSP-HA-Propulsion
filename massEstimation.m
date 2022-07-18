%% Structural Mass Estimate
% Estimates the mass of a Two-Stage rocket given geometires and masses for 
% airframes, fins and nosecons.  
%
% mass is reported in whatever units are desired as long as the units match 
% e.g density is lbs/in^3 and lengths are in inches. 

%% Fin Variables
clc
clear

%Fin_Material = 3; % Material of the fin

% FIRST STAGE
num_fins_1 = 4; % Number of fins on the first stage [in]
root_chord_1 = 14; % Length of fin where it meets the airframe [in]
tip_chord_1 = 5;% Length of the fin furthest from the airframe [in]
leading_angle_1 = 35; % Angle from the Root chord to leading edge [deg] (only used in CAD)
span_1 = 5; % Distance from Root Chord to Tip Chord [in]
fin_thickness_1 = .25; % Thickness of the fin [in]

% SECOND STAGE
num_fins_2 = 4; % Number of fins on the second stage [in]
root_chord_2 = 7; % Length of fin where it meets the airframe [in]
tip_chord_2 = 2.5; % Length of the fin furthest from the airframe [in]
leading_angle_2 = 35; % Angle from the Root chord to leading edge [deg] (only used in CAD)
span_2 = 3; % Distance from Root Chord to Tip Chord [in]
fin_thickness_2 = .25; % Thickness of the fin [in]

%% Airframe Variables

%AF_Material; % Material of the airframe

% FIRST STAGE
AF_id_1 = 3.9; % Inner Diameter of airframe [in]
AF_thickness_1 = 0.05; % Thickness of the airframe [in]
AF_height_1 = 45; % Hight of the airframe [in]

% SECOND STAGE
AF_id_2 = 3.9; % Inner Diameter of airframe [in]
AF_thickness_2 = 0.05; % Radial Thickness of the airframe [in] 
AF_height_2 = 30; % Hight of the airframe [in]

%% Nosecone Variables

%NC_Material; % Material of the hollow section of the nosecone
%Tip_Material; % Material of the tip of the nosecone

NC_shape = 2; %nosecone geometric shape [1-5]
NC_shape_parameter = 0; %modifies the nosecone geometry [relevant for haack series and power series];

NC_total_length = 20; % Overall length of the nosecone including the tip [in]
NC_thickness = 0.05; % Radial Thickness of the hollow section of the nosecone [in]
tip_length = 4.5; % Length of the solid tip of the nosecone [in]

%% Calculate Masses
% FINS
fin_density_1 = 0.1; % Mat Selection structure function [lb/in^3]
fin_density_2 = 0.1; % Mat Selection structure function [lb/in^3]

fin_mass_1 = num_fins_1 * finMass(fin_density_1, root_chord_1, tip_chord_1, span_1, fin_thickness_1); %first stage
fin_mass_2 = num_fins_2 * finMass(fin_density_2, root_chord_2, tip_chord_2, span_2, fin_thickness_2); %second stage
fin_mass = fin_mass_1 + fin_mass_2;

% AIRFRAMES
AF_Density = 0.07; % Mat Selection structure function [lb/in^3]

AF_mass_1 = airframeMass(AF_Density, AF_id_1, AF_thickness_1, AF_height_1);
AF_mass_2 = airframeMass(AF_Density, AF_id_2, AF_thickness_2, AF_height_2);
AF_mass = AF_mass_1 + AF_mass_2;

% NOSECONE

NC_density = 0.0556; % Mat Selection structure function [lb/in^3]
tip_density = 0.162; % Mat Selection structure function [lb/in^3]

NC_mass = noseconeMass(NC_shape, NC_shape_parameter, NC_thickness, NC_total_length, AF_thickness_2, AF_id_2, tip_density, NC_density, tip_length);

%TOTAL
structural_mass = fin_mass + AF_mass + NC_mass;
fprintf("Total esitmated mass is %d\n" ,structural_mass)

%% Export Variables to SolidWorks
