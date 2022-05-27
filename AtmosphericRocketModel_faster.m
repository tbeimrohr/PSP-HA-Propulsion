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
tic
clc
go = input('Run Trajectory Simulation? y = 1 n = 0\n');
if go
    clear
    go = 1;
    j = 0;
end

%% User Input Parameters
dv_change = 5.2384; %km/s
desired_alt = 300; %km
diameter1 = 4.5;%in
diameter2 = 4; %in
mpl = 1; %mass of the payload in kg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Model
dif_alt = 2;
alt_tol = .01;
check = 1;
o = 1;
while (abs(dif_alt - 1) > alt_tol) & check & o < 4 & go
    %% Constants
    j = j + 1;
    dt = .01; %sec
    At_srb = ((1/12)^2)*pi/39.37; %m2, the value is converted from inches^2 to m^2
    ep = [15]; %[26.86 29.2 22.1]; %[Formula#1 Formula#2 Formula#3 ...]
    At_srb2 = ((1/12)^2)*pi/39.37;  %m2, the value is converted from inches^2 to m^2
    ep2 = [10]; %[26.86 29.2 22.1]./3;

    pc_graph = [1000 975 980 990 995 985  925 875 810 785 735 690 680 655 650 625 625 650 650 400 375 250 10].*1.5; %chamber pressure history in psi
    %     pc_graph = [849	849.28125	849.28125	849.515625	849.796875	854.4375	854.625	855	838.3097015	814.9962687	784.2649254	784	784	688.3356643	634.6293706	564.1398601	544	544	401.836858	29.58308157	14.79154079	13.96978852	13.14803625	0]; %chamber pressure history in psi
    pc_graph2 = [600 650 675 700 710 715 720 725 730 735 740 745 750 750 750 750 750 750 750 750 750 600 100 10].*1.25; %chamber pressure history in psi

    stageHeight = [15]; %km, height at which staging must occur, this value is used in the interstage coast which can be turned off in the next line
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
        fprintf('Iteration %d...\n',j)
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

%% plots
figure(1)
    plot(t{j}(1:index_coast2(end)),F_srb{j}(1:index_coast2(end))./1000)
    hold on
    plot(t{j},Drag{j}./1000)
    plot(t{j},wDrag{j}./1000)
    plot(t{j},F_net{j}./1000)
    xlabel('Time [sec]')
    ylabel('Force [kN]')
    grid on
title('Forces Experienced up to 100 km')
legend('Rocket Thrust','Drag','Weight Drag','Net Force','location','best')



figure(2)
Details = [max(Q{j})/1000 t{j}(Q{j} == max(Q{j}))];
T2 = array2table(Details,'VariableNames',{'Maximum Dynamic Pressure [kPa]','Time of Max Q [sec]'});
T2 = table(T2,'VariableNames',{'Dynamic Pressure for Each Propellant'}); % Nested table
TString = evalc('disp(T2)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
set(gcf,'Position',[100 100 700 150])
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

figure(3)
plot(t{j},height{j}./1000)
hold on
plot(t{j},ones(1,length(t{j})).*desired_alt,'k--')
xlabel('Time [sec]')
ylabel('Altitude [km]')
grid on
title('Altitude of Rocket over Time')
%
% figure(4)
% subplot(1,2,1)
% plot(pc_graph)
% xlabel('Index Value')
% ylabel('Chamber Pressure [psi]')
% grid on
% title('Stage 1')
% subplot(1,2,2)
% plot(pc_graph2)
% xlabel('Index Value')
% ylabel('Chamber Pressure [psi]')
% grid on
% title('Stage 2')
% sgtitle('Chamber Pressure Profiles')

figure(5)
Details2 = [(mp{j}(1)+mp2{j}(1)) (mi(end)+mi2(end)) mo(end)];
T3 = array2table(Details2,'VariableNames',{'Propellant Mass [kg]','Inert Mass [kg]','Liftoff Mass [kg]'});
T3 = table(T3,'VariableNames',{'Total Launch Vehicle Masses Each Propellant'}); % Nested table
TString = evalc('disp(T3)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
set(gcf,'Position',[100 100 700 150])
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% figure(6)
% for i = 1:length(ve)
%     plot(t{j},acceleration{j}./g)
%     xlabel('Time [sec]')
%     ylabel('Acceleration [Gs]')
%     hold on
%     leg5{j} = strcat('Propellant #',string(i));
%     grid on
% end
% title('Acceleration of Rocket During Flight')
% legend(leg5,'location','best')

% figure(7)
% for i = 1:length(ve)
%     plot(t{j}(2:end),Mach{j})
%     xlabel('Time [sec]')
%     ylabel('Mach Number')
%     hold on
%     leg6{j} = strcat('Propellant #',string(i));
%     grid on
% end
% title('Mach Number of the Rocket During Flight')
% legend(leg6,'location','best')

figure(8)
Details3 = [tb(end) tb2(end) tot_dv_mission(end)];
T4 = array2table(Details3,'VariableNames',{'First Stage Burn Time [s]','Second Stage Burn Time [s]','dV [km/s]'});

T4 = table(T4,'VariableNames',{'Burn Time for Each Propellant'}); % Nested table
TString = evalc('disp(T4)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
set(gcf,'Position',[100 100 700 150])
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

toc