function [isp, cstar, ve, exit_pres, Tc] = runCEA(pc,ep)

addpath('C:\Users\1167948\Desktop\PSP GitHub\PSP-HA-Propulsion-CEA','-end');
savepath();

%% Input Weight Percentages of Formula
Al_wt = .15;
HTPB_wt = .1348;
AP_wt = .65;

%% Main Body of Code
%% CEA
r = 3.077;
g = 9.81;

inp = containers.Map;
inp('type') = 'eq';              % Sets the type of CEA calculation
inp('p') = pc;                % Chamber pressure
inp('p_unit') = 'psi';              % Chamber pressure units
inp('sup') = ep;
inp('fuel') = ["AL" "HTPB10"];             % Fuel name from thermo.inp
inp('fuel_wt%') = [Al_wt HTPB_wt];
inp('fuel_t') = [298 298];                % Fuel inlet temperature
inp('ox') = "NH4CLO4(I)";              % Ox name from thermo.inp
inp('ox_wt%') = AP_wt;
inp('ox_t') = 298;                  % Ox inlet temperature
inp('file_name') = sprintf('test1.inp');    % Input/output file name
inp('o/f') = r;               % Mixture ratio

data = cea_rocket_run(inp);     % Call the CEA MATLAB code

data_eq = data('eq');
% gamma = flip(squeeze(data_eq('gammas')));
% molmass = flip(squeeze(data_eq('m')));
Tc = flip(squeeze(data_eq('t')));
cstar = flip(squeeze(data_eq('cstar')));
cstar = cstar(:,1);
% eps = flip(squeeze(data_eq('ae/at')));
pressure = flip(squeeze(data_eq('p')));
exit_pres = pressure(:,end);
isps = flip(squeeze(data_eq('isp')));
isp = isps(:,end)./g;
mach = flip(squeeze(data_eq('mach')));
son = flip(squeeze(data_eq('son')));
ve = mach(:,end).*son(:,end);
