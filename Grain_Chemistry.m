function [data_srb] = Grain_Chemistry(pc)
% Change next line to match your folders
addpath('C:\Users\Thomas\Desktop\cea','-end');
savepath();

CEA_RUN = true;
CEA_SAVE_FILE = 'Project1.mat';

%% SRB
%% Input Weight Percentages of Formula
Al_wt = [.16];
PBAN_wt = [.14];
AP_wt = [.696];
FeO2_wt = [.004];


%% Main Body of Code
%% CEA
g = 9.81;
PBAN_molw = 99.9789;
AP_molw = 117.49;
Al_molw = 26.982;
FeO2_molw = 159.68820;

syms a b c d
systemEQ = [((a*PBAN_molw)/(a*PBAN_molw+b*AP_molw+c*Al_molw+d*FeO2_molw)) == PBAN_wt, ...
    ((b*AP_molw)/(a*PBAN_molw+b*AP_molw+c*Al_molw+d*FeO2_molw)) == AP_wt, ...
    ((c*Al_molw)/(a*PBAN_molw+b*AP_molw+c*Al_molw+d*FeO2_molw)) == Al_wt, ...
    ((d*FeO2_molw)/(a*PBAN_molw+b*AP_molw+c*Al_molw+d*FeO2_molw)) == FeO2_wt, ...
    a+b+c+d == 1];
mols = vpasolve(systemEQ,[a,b,c,d]);
A = double(mols.a);
B = double(mols.b);
C = double(mols.c);
D = double(mols.d);

mfuel = A.*PBAN_molw + C.*Al_molw;
mox = B.*AP_molw + D.*FeO2_molw;
r_srb = mox./mfuel;

re_srb = (152.6/2)/39.37; %m
rt_srb = (53.86/2)/39.37; %m
Ae_srb = pi*re_srb^2; %m2
At_srb = pi*rt_srb^2; %m2
ep_srb = Ae_srb/At_srb; %m2

for i = 1:length(Al_wt)
    inp = containers.Map;
    inp('type') = 'eq';              % Sets the type of CEA calculation
    inp('p') = pc;                % Chamber pressure
    inp('p_unit') = 'psi';              % Chamber pressure units
    inp('sup') = ep_srb;
    %inp('pip') = 1000/14.7;                   % Pressure ratios
    inp('fuel') = ["AL(cr)" "pban" "Fe2O3(S)"];             % Fuel name from thermo.inp
    inp('fuel_wt%') = [Al_wt(i) PBAN_wt(i) FeO2_wt(i)];
    inp('fuel_t') = [298 298 298];                % Fuel inlet temperature
    inp('ox') = "NH4CLO4(I)";              % Ox name from thermo.inp
    inp('ox_wt%') = AP_wt(i);
    inp('ox_t') = 298;                  % Ox inlet temperature
    inp('o/f') = r_srb(i);               % Mixture ratio
    inp('file_name') = sprintf('Project1.inp');    % Input/output file name
    if CEA_RUN
        data = cea_rocket_run(inp);     % Call the CEA MATLAB code
        save(CEA_SAVE_FILE, 'data');
    else
        load(CEA_SAVE_FILE);
    end
    data_srb = data('eq');

    gamma_srb = squeeze(data_srb('gammas'));

    molmass_srb = squeeze(data_srb('m'));

    Tc_srb = squeeze(data_srb('t'));

    cstar_srb = squeeze(data_srb('cstar'));

    eps_srb = squeeze(data_srb('ae/at'));

    pressure_srb = squeeze(data_srb('p'));

    isps_srb = squeeze(data_srb('isp'));
    isp_srb = isps_srb./g;

    mach_srb = squeeze(data_srb('mach'));

    son_srb = squeeze(data_srb('son'));
end
