clc
tic
addpath('C:\Users\Thomas\Desktop\cea','-end');
savepath();

CEA_RUN = true;
CEA_SAVE_FILE = 'TestCase.mat';

go = input("Run CEA? y = 1 n = 0\n");
if go
    clear
    go = 1;
    CEA_RUN = true;
    CEA_SAVE_FILE = 'test.mat';
    addpath('C:\Users\Thomas\Desktop\cea','-end');
    savepath();
    dw = 20;
    dw2 = 10;
    HTPB_molw = 99.99870;
    AP_molw = 117.49;
    Al_molw = 26.982;
    Al_wt = linspace(.1,.35,dw);
    HTPB_wt = linspace(.08,.2,dw2);

    for j = 1:length(HTPB_wt)
        for i = 1:length(Al_wt)
            AP_wt(j,i) = 1 - Al_wt(i) - HTPB_wt(j);
            syms a b c
            systemEQ = [((a*HTPB_molw)/(a*HTPB_molw+b*AP_molw+c*Al_molw)) == HTPB_wt(j), ...
                ((b*AP_molw)/(a*HTPB_molw+b*AP_molw+c*Al_molw)) == AP_wt(j,i), ...
                ((c*Al_molw)/(a*HTPB_molw+b*AP_molw+c*Al_molw)) == Al_wt(i), ...
                a+b+c == 1];
            mols = vpasolve(systemEQ,[a,b,c]);
            A = double(mols.a);
            B = double(mols.b);
            C = double(mols.c);

            mfuel = A.*HTPB_molw + C.*Al_molw;
            mox = B.*AP_molw;
            r(j,i) = mox./mfuel;
        end
    end
    g = 9.81;
    count = 0;
    HTPB_rho = 930;
    AP_rho = 1950;
    ep = 7.6;
end
if go
    for j = 1:length(HTPB_wt)
        for i = 1:length(Al_wt)
            count = count + 1;
            fprintf("Run %d of %d\n",count,length(Al_wt)*length(HTPB_wt))
            inp = containers.Map;
            inp('type') = 'eq';              % Sets the type of CEA calculation
            inp('p') = 1000;                % Chamber pressure
            inp('p_unit') = 'psi';              % Chamber pressure units
            inp('sup') = ep;               % Supersonic area ratios
            %inp('pip') = 1000/14.7;                   % Pressure ratios
            inp('fuel') = ["AL" "HTPB10"];             % Fuel name from thermo.inp
            inp('fuel_wt%') = [Al_wt(i) HTPB_wt(j)];
            inp('fuel_t') = [298 298];                % Fuel inlet temperature
            inp('ox') = "NH4CLO4(I)";              % Ox name from thermo.inp
            inp('ox_wt%') = AP_wt(j,i);
            inp('ox_t') = 298;                  % Ox inlet temperature
            inp('file_name') = sprintf('Design1%d%d.inp',i,j);    % Input/output file name
            inp('o/f') = r(j,i);               % Mixture ratio
            if CEA_RUN
                data = cea_rocket_run(inp);     % Call the CEA MATLAB code
                save(CEA_SAVE_FILE, 'data');
            else
                load(CEA_SAVE_FILE);
            end


            data_eq = data('eq');

            gamma = squeeze(data_eq('gammas'));
            Gammas(j,i) = gamma(1);

            molmass = squeeze(data_eq('m'));
            Molmass(j,i) = molmass(end);

            t = squeeze(data_eq('t'));
            Tc(j,i) = t(1);

            cstar = squeeze(data_eq('cstar'));
            Cstar(j,i) = cstar(1);

            isp = squeeze(data_eq('isp'));
            Isp(j,i) = isp(end)./g;
        end
    end
end
%% Plots
figure(1)
subplot(2,2,1)
sgtitle('CEA Investigation of HTPB, AP, & AL Composite Fuel')
plot(Al_wt,Molmass)
xlabel('Al Weight%')
ylabel('Molar Mass g/mol')
grid on
subplot(2,2,2)
plot(Al_wt,Gammas)
xlabel('Al Weight%')
ylabel('Gamma')
grid on
subplot(2,2,3)
plot(Al_wt,Tc)
xlabel('Al Weight%')
ylabel('Chamber Temp. K')
grid on
subplot(2,2,4)
plot(Al_wt,Cstar)
xlabel('Al Weight%')
ylabel('C^* m/s')
name = strcat(string(round(100.*HTPB_wt,2)),'% HTPB');

% add a bit space to the figure
set(gcf,'Position',[100 100 1500 500])
% add legend
Lgnd = legend(name);
Lgnd.Position(1) = .01;
Lgnd.Position(2) = 0.6;
grid on

figure(2)
plot(Al_wt,Isp)
xlabel('Al Weight%')
ylabel('Isp [sec]')
title('CEA Investigation of HTPB, AP, & AL Composite Fuel')
name = strcat(string(round(100.*HTPB_wt,2)),'% HTPB');
legend(name,'location','best')
grid on

toc