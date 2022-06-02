%% Atmospheric Rocket Model
%% Created by Thomas Beimrohr on 4/3/22
% For the purpose of creating a more defined model towards the development
% of a two stage rocket for PSP High Altitude.
%% Inputs
% Browse the section "Constants", in this section there are multiple
% constants to change. If there is a vector quantity, this means you can
% add more entries in order to model different propellents (use
% "Performance" or "Formulation" codes to obtain the required
% charecteristics of the propellant desired.
%% Start of Code
% tic
% clc
% go = input('Run Trajectory Simulation? y = 1 n = 0\n');
% if go
%     clear
%     close all

% end
function [score] = Atmospheric_1DoF(cpp1,cpp2,iter)
%% User Input Parameters
go = 1;
j = 0;
w = 0;
output_vars = {'Combination','Propellant Mass [kg]','Inert Mass [kg]',...
    'Total Mass [kg] (including payload)','First Stage Burn Time [sec]',...
    'Second Stage Burn Time [sec]','Maximum Dynamic Pressure (Max Q) [kPa]',...
    'Delta V [km/s]','Payload Mass (kg)','Alt [km]','Diameter 1','Diameter 2','Value','Cost','Score'};


out_table = array2table(zeros(1,length(output_vars)), 'VariableNames',output_vars);
row = 1;
col1 = 1;
col2 = 15;
RangeVariable1 = xlsAddr(row,col1);
RangeVariable2 = xlsAddr(row,col2);
RangeVariable = [RangeVariable1,':',RangeVariable2];
writetable(out_table,'ARM_Metrics.xls','Range',RangeVariable);
import_combo = readmatrix("Design Matrix - Praeto.csv",'Range', 'A2:A127','OutputType','string');
import_dv = readmatrix("Design Matrix - Praeto.csv", 'Range', 'H2:H127');

% for iter = (1:length(import_combo))
combo = import_combo(iter);
dv_change = import_dv(iter);

mpl_ref = combo{1}(1:2); %mass of the payload in kg
dia1_ref = combo{1}(3:4); %in
dia2_ref = combo{1}(5:6); %in
altitude_ref = combo{1}(7:8); %km


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% CSV Identification Codes
mpl_code = {{'01' '02' '03'} [1 5 10]};
dia1_code = {{'04' '05' '06'} [4 4.5 5]};
dia2_code = {{'07' '08'}, [3.5 4]};
altitude_code = {{'09' '10' '11' '12' '13' '14' '15'}, [100 125 150 200 250 300 400]};

mpl = mpl_code{2}(find(strcmp(mpl_ref, mpl_code{1}) == 1));
diameter1 = dia1_code{2}(find(strcmp(dia1_ref, dia1_code{1}) == 1));
diameter2 = dia2_code{2}(find(strcmp(dia2_ref, dia2_code{1}) == 1));
desired_alt = altitude_code{2}(find(strcmp(altitude_ref, altitude_code{1}) == 1));


%% Model
dif_alt = 2;
alt_tol = .01;
check = 1;
o = 1;
while (abs(dif_alt - 1) > alt_tol) & check & o < 4 & go
    %% Constants
    j = j + 1;
    w = w + 1;
    dt = .01; %sec
    At_srb = ((1/12)^2)*pi/39.37; %m2, the value is converted from inches^2 to m^2
    ep = [15]; %[26.86 29.2 22.1]; %[Formula#1 Formula#2 Formula#3 ...]
    At_srb2 = ((1/12)^2)*pi/39.37;  %m2, the value is converted from inches^2 to m^2
    ep2 = [10]; %[26.86 29.2 22.1]./3;

    pc_graph = cpp1.coords;
    pc_graph2 = cpp2.coords;

    stageHeight = [22]; %km, height at which staging must occur, this value is used in the interstage coast which can be turned off in the next line
    BeforeStageCoast = 0; %value used to turn off or on the interstage coast, bewarned if the value of stageHeight is chosen to be too high where the rocket doesnt make it to the desired height, the answers obtained will become nonsensical

    cstar_srb = [1737.5]; %m/s
    ve = [2461.1]; %m/s
    Pe_srb = [29073]; %Pa
    Pe_srb2 = [29073]; %Pa

    tol = .00001; %this is the tolerance used to calculate propellant mass vs the given propellant mass, the code converges rather quickly so this value can be set extremely small to obtain "more accurate" results
    difmp = ones(1,length(ve)).*10; %initalizing the calculated prop mass difference
    tb(j) = ones(1,length(ve)).*1; %guessing the burn out time of stage 1
    difmp2 = ones(1,length(ve)).*10; %initalizing the calculated prop mass difference
    tb2(j) = ones(1,length(ve)).*1; %guessing the burn out time of stage 2

    lam1 = .85; %propellant mass fraction of stage 1
    lam2 = .785; %propellant mass fraction of stage 2


    Isp = [280].*.9; %sec
    g = 9.81; %m/s2
    tot_dv_mission(j) = [dv_change]; %km/s
    scale1 = .37; %the amount of delta v percentage the first stage will carry
    Area1 = (pi*(diameter1/2)^2)/1550; %m2
    Area2 = (pi*(diameter2/2)^2)/1550; %m2

    Cd_o = .5;
    %% 1D Model
    Ae_srb(j) = ep * At_srb;
    Ae_srb2(j) = ep2 * At_srb2;

    dv1 = scale1.*tot_dv_mission(j);
    dv2 = (1-scale1).*tot_dv_mission(j);

    MR1 = exp((dv1.*1000)./(g.*Isp));
    MR2 = exp((dv2.*1000)./(g.*Isp));

    mp2{j} = mpl.*((MR2-1)./(MR2-((MR2-1)./lam2))); %kg
    mi2(j) = mp2{j}.*((1-lam2)./lam2); %kg
    mo2(j) = mi2(j) + mpl + mp2{j}; %kg
    mf2(j) = mi2(j) + mpl; %kg

    mp{j} = mo2(j).*((MR1-1)./(MR1-((MR1-1)./lam1))); %kg
    mi(j) = mp{j}.*((1-lam1)./lam1); %kg
    mo(j) = mi(j) + mo2(j) + mp{j}; %kg
    mf(j) = mi(j) + mo2(j); %kg

    while abs(difmp - 1) > tol
        time_srb = linspace(0,tb(j),length(pc_graph)); %sec
        pc_srb = ModifySize(pc_graph,length(time_srb)); %psi
        Pc_srb = pc_srb.*6894.76; %Pa
        mdot_srbA{j} = Pc_srb.*At_srb./cstar_srb;
        intmdot = cumtrapz(mdot_srbA{j})*(time_srb(2)-time_srb(1));
        mdot_check = intmdot(end);

        difmp = mdot_check/mp{j};
        tb(j) = tb(j) + 2*(1 - difmp);
    end

    time{j} = linspace(0,time_srb(end),time_srb(end)/dt + 1);
    mdot_srb{j} = interp1(time_srb,mdot_srbA{j},time{j});
    mdotv_srb{j} = mdot_srb{j}.*ve;

    while abs(difmp2 - 1) > tol
        time_srb2 = linspace(0,tb2(j),length(pc_graph2)); %sec
        pc_srb2 = ModifySize(pc_graph2,length(time_srb2)); %psi
        Pc_srb2 = pc_srb2.*6894.76; %Pa
        mdot_srbB{j} = Pc_srb2.*At_srb2./cstar_srb;
        intmdot2 = cumtrapz(mdot_srbB{j})*(time_srb2(2)-time_srb2(1));
        mdot_check2 = intmdot2(end);

        difmp2 = mdot_check2/mp2{j};
        tb2(j) = tb2(j) + 1*(1 - difmp2);
    end

    time2{j} = linspace(0,time_srb2(end),time_srb2(end)/dt + 1);
    mdot_srb2{j} = interp1(time_srb2,mdot_srbB{j},time2{j});
    mdotv_srb2{j} = mdot_srb2{j}.*ve;

    i = 1;
    t{j}(1) = 0; %sec
    F_srb{j}(1) = 0;
    vel{j}(1) = 0; %m/s
    m{j}(1) = mo(j);
    height{j}(1) = 0; %*m
    [T{j}(1),Son{j}(1),P{j}(1),Rho{j}(1)] = atmoscoesa(height{j}(1),'None');

    if go
        clc
        fprintf('Combination %d...\n',iter)
        fprintf('Iteration %d...\n',w)
        fprintf('First Stage Burn...\n')
        for i = 1:length(time{j})
            t{j}(i+1) = t{j}(i) + dt;
            [T{j}(i),Son{j}(i),P{j}(i),Rho{j}(i)] = atmoscoesa(height{j}(i),'None');
            Mach{j}(i) = vel{j}(i)/Son{j}(i);
            mp{j}(i+1) = mp{j}(i) - mdot_srb{j}(i)*dt;
            m{j}(i+1) = mp{j}(i+1) + mi(j) + mo2(j);
            F_srb{j}(i+1) = (mdotv_srb{j}(i) + Ae_srb(j).*(Pe_srb - P{j}(i)));
            if Mach{j}(i) < .5
                Cd{j}(i) = Cd_o;
                %                 Cd{j}(i) = Cd_o/(sqrt(1-(Mach{j}(i))^2));
            elseif Mach{j}(i) >= .5 & Mach{j}(i) <= 2 & Cd{j}(i-1) > Cd_o*.97
                Cd{j}(i) = Cd_o^((1-Mach{j}(i))*sign(1-Mach{j}(i)));
            else
                Cd{j}(i) = Cd_o*.97;
            end
            Drag{j}(i+1) = .5*Rho{j}(i)*Cd{j}(i)*(vel{j}(i)^2)*Area1;
            wDrag{j}(i+1) = m{j}(i)*g;
            F_net{j}(i+1) = F_srb{j}(i) - Drag{j}(i) - wDrag{j}(i);
            acceleration{j}(i+1) = F_net{j}(i)/m{j}(i);
            if i >= 1
                intacc = vel{j}(i) + cumtrapz(acceleration{j}(end-1:end))*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = height{j}(i) + cumtrapz(vel{j}(end-1:end))*dt;
                height{j}(i+1) = intvel(end);
            else
                intacc =cumtrapz(acceleration{j})*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = cumtrapz(vel{j})*dt;
                height{j}(i+1) = intvel(end);
            end
            i = i + 1;

        end
        index_coast(j) = i;
        fprintf('Interstage Coast...\n')
        while height{j}(i)/1000 <= stageHeight & BeforeStageCoast
            t{j}(i+1) = t{j}(i) + dt;
            [T{j}(i),Son{j}(i),P{j}(i),Rho{j}(i)] = atmoscoesa(height{j}(i),'None');
            if isnan(P{j}(i))
                P{j}(i) = 0;
                Rho{j}(i) = 0;
            end
            Mach{j}(i) = vel{j}(i)/Son{j}(i);
            m{j}(i+1) = m{j}(end);
            if Mach{j}(i) < .5
                Cd{j}(i) = Cd_o;
                %                 Cd{j}(i) = Cd_o/(sqrt(1-(Mach{j}(i))^2));
            elseif Mach{j}(i) >= .5 & Mach{j}(i) <= 2 & Cd{j}(i-1) > Cd_o*.97
                Cd{j}(i) = Cd_o^((1-Mach{j}(i))*sign(1-Mach{j}(i)));
            else
                Cd{j}(i) = Cd_o*.97;
            end
            Drag{j}(i+1) = .5*Rho{j}(i)*Cd{j}(i)*(vel{j}(i)^2)*Area1;
            wDrag{j}(i+1) = m{j}(i)*g;
            F_net{j}(i+1) = -wDrag{j}(i) - Drag{j}(i);
            acceleration{j}(i+1) = F_net{j}(i+1)/m{j}(end);
            if i >= 1+index_coast(j)
                intacc = vel{j}(i) + cumtrapz(acceleration{j}(end-1:end))*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = height{j}(i) + cumtrapz(vel{j}(end-1:end))*dt;
                height{j}(i+1) = intvel(end);
            else
                intacc =cumtrapz(acceleration{j})*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = cumtrapz(vel{j})*dt;
                height{j}(i+1) = intvel(end);
            end
            i = i + 1;
        end
        index_stage(j) = i;
        k = 1;
        fprintf('Second Stage Stage Burn...\n')
        for i = index_stage(j):index_stage(j) + length(time2{j})-1
            t{j}(i+1) = t{j}(i) + dt;
            [T{j}(i),Son{j}(i),P{j}(i),Rho{j}(i)] = atmoscoesa(height{j}(i),'None');
            if isnan(P{j}(i))
                P{j}(i) = 0;
                Rho{j}(i) = 0;
            end
            Mach{j}(i) = vel{j}(i)/Son{j}(i);
            if Mach{j}(i) < .5
                Cd{j}(i) = Cd_o;
                %                 Cd{j}(i) = Cd_o/(sqrt(1-(Mach{j}(i))^2));
            elseif Mach{j}(i) >= .5 & Mach{j}(i) <= 2 & Cd{j}(i-1) > Cd_o*.97
                Cd{j}(i) = Cd_o^((1-Mach{j}(i))*sign(1-Mach{j}(i)));
            else
                Cd{j}(i) = Cd_o*.97;
            end
            mp2{j}(k+1) = mp2{j}(k) - mdot_srb2{j}(k)*dt;
            m{j}(i+1) = mp2{j}(k+1) + mi2(j) + mpl;
            F_srb{j}(i+1) = (mdotv_srb2{j}(k) + Ae_srb2(j).*(Pe_srb2 - P{j}(i)));
            Drag{j}(i+1) = .5*Rho{j}(i)*Cd{j}(i)*(vel{j}(i)^2)*Area2;
            wDrag{j}(i+1) = m{j}(i)*g;
            F_net{j}(i+1) = F_srb{j}(i) - Drag{j}(i) - wDrag{j}(i);
            acceleration{j}(i+1) = F_net{j}(i)/m{j}(i);
            if i >= 1+index_stage(j)
                intacc = vel{j}(i) + cumtrapz(acceleration{j}(end-1:end))*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = height{j}(i) + cumtrapz(vel{j}(end-1:end))*dt;
                height{j}(i+1) = intvel(end);
            else
                intacc =cumtrapz(acceleration{j})*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = cumtrapz(vel{j})*dt;
                height{j}(i+1) = intvel(end);
            end
            i = i + 1;
            k = k + 1;
        end
        index_coast2(j) = i;
        fprintf('Final Coast...\n')
        while vel{j}(i) >= 0
            t{j}(i+1) = t{j}(i) + dt;
            [T{j}(i),Son{j}(i),P{j}(i),Rho{j}(i)] = atmoscoesa(height{j}(i),'None');
            if isnan(P{j}(i))
                P{j}(i) = 0;
                Rho{j}(i) = 0;
            end
            Mach{j}(i) = vel{j}(i)/Son{j}(i);
            if Mach{j}(i) < .5
                Cd{j}(i) = Cd_o;
                %                 Cd{j}(i) = Cd_o/(sqrt(1-(Mach{j}(i))^2));
            elseif Mach{j}(i) >= .5 & Mach{j}(i) <= 2 & Cd{j}(i-1) > Cd_o*.97
                Cd{j}(i) = Cd_o^((1-Mach{j}(i))*sign(1-Mach{j}(i)));
            else
                Cd{j}(i) = Cd_o*.97;
            end
            m{j}(i+1) = m{j}(end);
            Drag{j}(i+1) = .5*Rho{j}(i)*Cd{j}(i)*(vel{j}(i)^2)*Area2;
            wDrag{j}(i+1) = m{j}(i)*g;
            F_net{j}(i+1) = -wDrag{j}(i) - Drag{j}(i);
            acceleration{j}(i+1) = F_net{j}(i+1)/m{j}(end);
            if i >= 1+index_coast2(j)
                intacc = vel{j}(i) + trapz(acceleration{j}(end-1:end)).*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = height{j}(i) + trapz(vel{j}(end-1:end)).*dt;
                height{j}(i+1) = intvel(end);
            else
                intacc = cumtrapz(acceleration{j})*dt;
                vel{j}(i+1) = intacc(end);
                Q{j}(i+1) = .5*Rho{j}(i)*vel{j}(i)^2;
                intvel = cumtrapz(vel{j})*dt;
                height{j}(i+1) = intvel(end);
            end
            i = i + 1;
        end
    end

    dif_alt = (height{j}(end)/(desired_alt*1000));
    sign_check(j) = sign(dif_alt - 1);
    if j > 1 && sign_check(j) ~= sign_check(j-1)
        o = o + 1;
    end
    dv_change = dv_change + (.1/(o^2))*((1-dif_alt)/dif_alt);
end
%% Ouputs
identif = str2double(append(mpl_ref, dia1_ref, dia2_ref, altitude_ref));
output_vec = {identif, (mp{j}(1)+mp2{j}(1)),(mi(end)+mi2(end)), mo(end), ...
    tb(end), tb2(end), max(Q{j})/1000, tot_dv_mission(end), mpl, desired_alt, diameter1, diameter2};

[~,Value,Cost] = Pareto_Selection(output_vec);
score = Value - Cost;

output_vec = {identif, (mp{j}(1)+mp2{j}(1)),(mi(end)+mi2(end)), mo(end), ...
    tb(end), tb2(end), max(Q{j})/1000, tot_dv_mission(end), mpl, desired_alt, diameter1, diameter2,Value,Cost, score};

out_table = array2table(output_vec);
row = iter+1;
col1 = 1;
col2 = 15;
RangeVariable1 = xlsAddr(row,col1);
RangeVariable2 = xlsAddr(row,col2);
RangeVariable = [RangeVariable1,':',RangeVariable2];
writematrix(cell2mat(output_vec),'ARM_Metrics.xls','Range',RangeVariable);


hold off
%     figure(2)
%     plot(t{j},height{j})
%     fprintf("%.6f",height{j}(end))
end
