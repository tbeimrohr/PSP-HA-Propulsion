clc
tic
% Change next line to match your folders
addpath('C:\Users\Thomas\Desktop\cea','-end');
savepath();

CEA_RUN = true;
CEA_SAVE_FILE = 'test.mat';
go = input("Run CEA? y = 1 n = 0\n");
if go
    clear
    accuracy = input("How accurate? 1 thru 10\n");
    go = 1;
    CEA_RUN = true;
    CEA_SAVE_FILE = 'test.mat';
    addpath('C:\Users\Thomas\Desktop\cea','-end');
    savepath();
end

%% Input Weight Percentages of Formula
Al_wt = [.25];
HTPB_wt = [.133333333];
AP_wt = 1-Al_wt-HTPB_wt;

%% Main Body of Code
%% CEA
HTPB_molw = 99.99870;
AP_molw = 117.49;
Al_molw = 26.982;

syms a b c
systemEQ = [((a*HTPB_molw)/(a*HTPB_molw+b*AP_molw+c*Al_molw)) == HTPB_wt, ...
    ((b*AP_molw)/(a*HTPB_molw+b*AP_molw+c*Al_molw)) == AP_wt, ...
    ((c*Al_molw)/(a*HTPB_molw+b*AP_molw+c*Al_molw)) == Al_wt];
mols = vpasolve(systemEQ,[a,b,c]);
A = double(mols.a);
B = double(mols.b);
C = double(mols.c);

mfuel = A.*HTPB_molw + C.*Al_molw;
mox = B.*AP_molw;
r = mox./mfuel;

g = 9.81;
HTPB_rho = 930;
AP_rho = 1950;
pc = 1000;
pa = [14.7 14.7];
ep(1,:) = linspace(0,7.6,6*accuracy);
for i = 1:length(Al_wt)
    if go
        fprintf('Run %d out of %d\n',i,length(Al_wt))
        inp = containers.Map;
        inp('type') = 'eq';              % Sets the type of CEA calculation
        inp('p') = pc;                % Chamber pressure
        inp('p_unit') = 'psi';              % Chamber pressure units
        inp('sup') = ep(i,:);
        %inp('pip') = 1000/14.7;                   % Pressure ratios
        inp('fuel') = ["AL" "HTPB10"];             % Fuel name from thermo.inp
        inp('fuel_wt%') = [Al_wt(i) HTPB_wt(i)];
        inp('fuel_t') = [298 298];                % Fuel inlet temperature
        inp('ox') = "NH4CLO4(I)";              % Ox name from thermo.inp
        inp('ox_wt%') = AP_wt(i);
        inp('ox_t') = 298;                  % Ox inlet temperature
        inp('file_name') = sprintf('test%d.inp',i);    % Input/output file name
        inp('o/f') = r(i);               % Mixture ratio
        if CEA_RUN
            data = cea_rocket_run(inp);     % Call the CEA MATLAB code
            save(CEA_SAVE_FILE, 'data');
        else
            load(CEA_SAVE_FILE);
        end
        data_eq = data('eq');
        
        gamma(:,i) = squeeze(data_eq('gammas'));
        
        molmass(:,i) = squeeze(data_eq('m'));
        
        Tc(:,i) = squeeze(data_eq('t'));
        
        cstar(:,i) = squeeze(data_eq('cstar'));
        
        eps(:,i) = squeeze(data_eq('ae/at'));
        
        pressure(:,i) = squeeze(data_eq('p'));
        
        isps(:,i) = squeeze(data_eq('isp'));
        isp = isps./g;

        mach(:,i) = squeeze(data_eq('m'));
        son(:,i) = squeeze(data_eq('son'));
    end
end
%% Calculations
ve = mach(end)*son(end)
pressure(end)
for i = 1:length(Al_wt)
    % nozzle seperation pressure
    Psep = pa(i)*.667*((pc/pa(1))^(-.2));
    Psep_pa = Psep*6894.76;
    Pexit = pressure(end,i);
    if Pexit < Psep_pa
        pressures = pressure(1,i):-1:pressure(end,i);
        vq = interp1(pressure(:,i),eps,pressures);
        dif = abs(pressures - Psep_pa);
        indexdif = find(dif == min(dif));
        eps_closest = vq(indexdif);
        sprintf('bad, seperation at exit for choice %d!\nepsilon of %f is about right',i,eps_closest)
    end
end

%% Plots
figure(1)
subplot(3,2,1)
sgtitle({'Flow Properties of 13.3% HTPB, 25% AL, 61.7% AP','at a Chamber Pressure of 1000 psi'})
plot(eps(:,1),molmass(:,1))
% hold on
% plot(eps(:,2),molmass(:,2))
xlabel('Expansion Ratio')
ylabel('Molar Mass [g/mol]')
grid on
% legend('New Formula','Old Formula','location','best')
subplot(3,2,2)
plot(eps(:,1),gamma(:,1))
% hold on
% % plot(eps(:,2),gamma(:,2))
xlabel('Expansion Ratio')
ylabel('Specific Heat Ratio \gamma')
grid on
subplot(3,2,3)
plot(eps(:,1),Tc(:,1))
% hold on
% plot(eps(:,2),Tc(:,2))
xlabel('Expansion Ratio')
ylabel('Temperature [K]')
grid on
% subplot(2,2,3)
% plot(eps(:,1),cstar(:,1))
% % hold on
% % plot(eps(:,2),cstar(:,2))
xlabel('Expansion Ratio')
% ylabel('C^* [m/s]')
% grid on
subplot(3,2,4)
plot(eps(:,1),isp(:,1))
% hold on
% plot(eps(:,2),isp(:,2))
xlabel('Expansion Ratio')
ylabel('Isp [sec]')
grid on
subplot(3,2,5)
plot(eps(:,1),(pressure(:,1)./6895))
% hold on
% plot(eps(:,2),pressure(:,2))
xlabel('Expansion Ratio')
ylabel('Pressure [psi]')
grid on

% figure(2)
% plot(eps(:,1),isp(:,1))
% hold on
% plot(eps(:,2),isp(:,2))
% xlabel('Expansion Ratio')
% ylabel('Isp sec')
% legend('New Formula','Old Formula','location','best')
% grid on
% 
% figure(3)
% Choice = [1;length(Al_wt)];
% Max_Isp = isp(end,:)';
% T = table(Choice,Max_Isp);
% TString = evalc('disp(T)');
% TString = strrep(TString,'<strong>','\bf');
% TString = strrep(TString,'</strong>','\rm');
% TString = strrep(TString,'_','\_');
% FixedWidth = get(0,'FixedWidthFontName');
% set(gcf,'Position',[300 300 250 100])
% annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
%     'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
toc