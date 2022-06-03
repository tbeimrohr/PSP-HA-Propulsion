%% Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:
% RunGeneticTrial
%
% Date of Creation:
% 04/28/2022
%
% Author(s):
% Thomas Beimrohr
% Jeff Kaji
% Justin Kruse
%
% Description:
% During the Pareto selection process, the need for an adaptable system
% that could test a set list of design parameter combinations was needed.
% This program is the master code that will run the 1 dimensional model and
% fill that model in with a genetic trial to produce optimal thrust
% profiles for the first and second stage of the rocket.
%
% Inputs:
% (required) pc_graph = the initial seed thrust profile for the first stage
% (required) pc_graph2 = the initial seed thrust profile for the second stage
% (required) import_combo = a .csv that contains information about the
% rocket stored in the appropriate formating standards
%
% Outputs:
% (required) ARM_Metrics.csv = .csv file with key information: Combination ID
% Propellant Mass [kg]	Inert Mass [kg]	Total Mass [kg] (including payload)	
% First Stage Burn Time [sec]	Second Stage Burn Time [sec]	
% Maximum Dynamic Pressure (Max Q) [kPa]	Delta V [km/s]	
% Payload Mass (kg)	Alt [km]	Diameter 1	Diameter 2
%
%
% Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of Code
clear
clc
tic
pc_graph = [1000 975 980 990 995 985 925 875 810 785 735 690 680 655 650 625 625 650 650 400 375 250 10 0].*1.25; %chamber pressure history in psi
pc_graph2 = [600 650 675 700 710 715 720 725 730 735 740 745 750 750 750 750 750 750 750 750 600 100 10 0].*1.25; %chamber pressure history in psi
rng(5)
id = 1; % Combonation id sequencial order
numBest = 1; % Take best n number of profiles from each generation
numGens = 1; % Generations to run for

%import_combo = readmatrix("Design Matrix - Praeto.csv",'Range', 'A2:A127','OutputType','string');
Pareto_write
a = string(load('combos.txt'));
import_combo = replace(a, '-', '0');

for i = 1:length(import_combo)
    b = run(i,numBest,numGens,curve1=pc_graph,curve2=pc_graph2,stage=1);
end
toc
