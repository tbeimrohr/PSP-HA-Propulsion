%% Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:
% AtmosphericRocketModel_faster
%
% Date of Creation:
% 04/3/2022
%
% Author(s):
% Thomas Beimrohr
% Jeff Kaji
% Justin Kruse
%
% Description:
% For the purpose of creating a more defined model towards the development
% of a two stage rocket for PSP High Altitude. The model is a 1 dimensional
% model that takes into account Cd variation, net forces of thrust drag
% gravitational drag. The code as a whole is used to size the propulsion
% system and then using empirical data size the inert mass of the rocket
%
% Inputs:
% (required) dv_change = initial guess for the amount of delta v mission
% will take
% (required) desired_alt = desired altitude obtained by the rocket
% (required) diameter1 = diameter of the first stage of the rocket
% (required) diameter2 = diameter of the second stage of the rocket
% (required) mpl = mass of the payload
% (required) At_srb = area of the throat of the first stage of the rocket
% (required) ep = expansion ratio of first stage
% (required) At_srb2 = area of the throat of the second stage of the rocket
% (required) ep2 = expansion ratio of second stage
% (required) pc_graph = chamber pressure history for the first stage
% (required) pc_graph2 = chamber pressure history for the second stage
% (required) stageHeight = height at which seperation of first and second
% stage should occur
% (required) BeforeStageCoast = turns on or off the coasting function tied
% to stageHeight
% (required) cstar_srb = the c* for the propellant being used for both
% first and second stage
% (required) ve = exit velocity from the nozzle
% (required) Pe_srb = exit pressure from the first stage nozzle
% (required) Pe_srb2 = exit pressure from the first stage nozzle
% (required) lam1 = propellant mass fraction of first stage from empirical data
% lambda = (mass of propellant)/(total mass of system)
% (required) lam2 = propellant mass fraction of second stage from empirical data
% (required) Isp = Isp of the propellant chosen
% (required) scale1 = delta v spread between the first and second stages
% (required) Cd_o = coeffiencient of drag for rocket body
%
% Outputs:
% (optional) Drag = time history for the drag force
% (optional) Fnet = time history for the net force acting on the rocket
% (optional) F_srb = time history for the thrust of the rocket
% (optional) height = time history for the altitude obtained by the rocket
% (optional) Mach = time history for the Mach of the rocket
% (optional) mdot_srb = time history for the mass flowrate from stage 1
% (optional) mdot_srb2 = time history for the mass flow rate of stage 2
% (optional) m = time history for the total mass of the rocket
% (optional) P = time history of atmospheric pressure
% (optional) Q = time history for the dynamic pressure
% (optional) Rho = time history of the atmospheric density
% (optional) Son = time history of the atmospheric sonic velocity
% (optional) T = time history of the atmospheric Temperature
% (optional) t = time history for the time o flight
% (optional) tb1(end) = burnout time of the first stage
% (optional) tb2(end) = burnout time of the second stage
% (optional) vel = time history for the velocity of the rocket
% (optional) wDrag = time history for the weight drag of the rocket
% (optional) acceleration = time history for the acceleration of the rocket
%
% Notes:
% - Code is sensitive to the starting guess for delta v in "dv_change", if
% code doesnt run aka time = nan, then try changing dv_change. 
% - If the coasting option is turned on, and answers become either unreasonable or
% nan or doesnt converge then the issue is the coasting height set in
% "StageHeight" is too high, ie the rocket is stopped by drag and gravity
% before it can reach the set point for second stage ignition. 
% - constants with brackets [] are used for testing multiple propellants at
% the same time
% - The only known improvement at this time is to corelate the exit
% pressure of the nozzle to the chamber pressure profile added by
% "pc_graph" this can be done within CEA. Possibly add a CEA function to
% produce this value and or the rest of the propellant inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of Code
tic
clc
go = input('Run Trajectory Simulation? y = 1 n = 0\n');
if go
    clear
    go = 1;
    j = 0;
end
%% Constants
dv_change = 5.097517151; %km/s
desired_alt = 100; %km
diameter1 = 4.5;%in
diameter2 = 4.5; %in
mpl = 1; %mass of the payload in kg
scale1 = .35; %the amount of delta v percentage the first stage will carry
%% Model
dif_alt = 2;
alt_tol = .01;
check = 1;
o = 1;
if ~exist("CEA_PcData_15.mat")
    pc1_cea = 1:25:1750;
    [isp1_cea, cstar1_cea, ve1_cea, exit_pressure1_cea, Tc1_cea] = runCEA(pc1_cea,ep);
    save("CEA_PcData_15.mat","isp1_cea","cstar1_cea","ve1_cea","exit_pressure1_cea","Tc1_cea","pc1_cea")
else
    load("CEA_PcData_15");
end
if ~exist("CEA_PcData_10.mat")
    pc2_cea = 1:20:1100;
    [isp2_cea, cstar2_cea, ve2_cea, exit_pressure2_cea, Tc2_cea] = runCEA(pc2_cea,ep2);
    save("CEA_PcData_10.mat","isp2_cea","cstar2_cea","ve2_cea","exit_pressure2_cea","Tc2_cea","pc2_cea")
else
    load("CEA_PcData_10");
end
while (abs(dif_alt - 1) > alt_tol) & check & o < 4 & go
    %% Constants 2
    j = j + 1;
    dt = .01; %sec
    At_srb = ((1/12)^2)*pi/39.37; %m2, the value is converted from inches^2 to m^2
    ep = [15]; %[26.86 29.2 22.1]; %[Formula#1 Formula#2 Formula#3 ...]
    At_srb2 = ((1/12)^2)*pi/39.37;  %m2, the value is converted from inches^2 to m^2
    ep2 = [10]; %[26.86 29.2 22.1]./3;

    pc_graph = [1433	1433	1444	1451	1454	1468	1470	1460	1435	1332	1295	1139	1138	1065	720	675	642	637	607	506	480	472	415	0]; %chamber pressure history in psi
    pc_graph2 = [643	668.5	705.125	719.8125	722	722	722	721	721.2857143	720	680.5714286	720.2745098	602.4444444	712.3333333	788	797.1687243	726	731.308642	704.1481481	653.5714286	593.4151786	464.9732143	214.9616071	23.2]; %chamber pressure history in psi

    pc_graph(pc_graph < 1) = 1;
    pc_graph2(pc_graph2 < 1) = 1;

    isp1 = interp1(pc1_cea,isp1_cea,pc_graph);
    cstar1 = interp1(pc1_cea,cstar1_cea,pc_graph);
    ve1 = interp1(pc1_cea,ve1_cea,pc_graph);
    exit_pressure1 = interp1(pc1_cea,exit_pressure1_cea,pc_graph);
    Tc1 = interp1(pc1_cea,Tc1_cea,pc_graph);

    isp2 = interp1(pc2_cea,isp2_cea,pc_graph2);
    cstar2 = interp1(pc2_cea,cstar2_cea,pc_graph2);
    ve2 = interp1(pc2_cea,ve2_cea,pc_graph2);
    exit_pressure2 = interp1(pc2_cea,exit_pressure2_cea,pc_graph2);
    Tc2 = interp1(pc2_cea,Tc2_cea,pc_graph2);

    Pe_srb = exit_pressure1; %Pa
    Pe_srb2 = exit_pressure2; %Pa

    stageHeight = [15]; %km, height at which staging must occur, this value is used in the interstage coast which can be turned off in the next line
    BeforeStageCoast = 0; %value used to turn off or on the interstage coast, bewarned if the value of stageHeight is chosen to be too high where the rocket doesnt make it to the desired height, the answers obtained will become nonsensical

    tol = .00001; %this is the tolerance used to calculate propellant mass vs the given propellant mass, the code converges rather quickly so this value can be set extremely small to obtain "more accurate" results
    difmp = 10; %initalizing the calculated prop mass difference
    tb(j) = 1; %guessing the burn out time of stage 1
    difmp2 = 10; %initalizing the calculated prop mass difference
    tb2(j) = 1; %guessing the burn out time of stage 2

    lam1 = .85; %propellant mass fraction of stage 1
    lam2 = .785; %propellant mass fraction of stage 2


    eta_isp = .925;
    Isp1 = mean(isp1)*eta_isp; %sec
    Isp2 = mean(isp2)*eta_isp;
    g = 9.81; %m/s2
    tot_dv_mission(j) = [dv_change]; %km/s
    Area1 = (pi*(diameter1/2)^2)/1550; %m2
    Area2 = (pi*(diameter2/2)^2)/1550; %m2

    Cd_o = .5;
    %% 1D Model
    Ae_srb(j) = ep * At_srb;
    Ae_srb2(j) = ep2 * At_srb2;

    dv1 = scale1.*tot_dv_mission(j);
    dv2 = (1-scale1).*tot_dv_mission(j);

    MR1 = exp((dv1.*1000)./(g.*Isp1));
    MR2 = exp((dv2.*1000)./(g.*Isp2));

    mp2{j} = mpl.*((MR2-1)./(MR2-((MR2-1)./lam2))); %kg
    mi2(j) = mp2{j}.*((1-lam2)./lam2); %kg
    mo2(j) = mi2(j) + mpl + mp2{j}; %kg
    mf2(j) = mi2(j) + mpl; %kg

    mp{j} = mo2(j).*((MR1-1)./(MR1-((MR1-1)./lam1))); %kg
    mi(j) = mp{j}.*((1-lam1)./lam1); %kg
    mo(j) = mi(j) + mo2(j) + mp{j}; %kg
    mf(j) = mi(j) + mo2(j); %kg

    while abs(difmp - 1) > tol
        tb(j) = tb(j) + 2*(1-difmp);
        time_srb = linspace(0,tb(j),length(pc_graph)); %sec
        pc_srb = ModifySize(pc_graph,length(time_srb)); %psi
        cstar1_srb = ModifySize(cstar1,length(time_srb));
        Pc_srb = pc_srb.*6894.76; %Pa
        mdot_srbA{j} = Pc_srb.*At_srb./cstar1_srb;
        intmdot = cumtrapz(mdot_srbA{j})*(time_srb(2)-time_srb(1));
        mdot_check = intmdot(end);

        difmp = mdot_check/mp{j};
    end

    time{j} = linspace(0,time_srb(end),time_srb(end)/dt + 1);
    mdot_srb{j} = interp1(time_srb,mdot_srbA{j},time{j});
    ve1_srb = ModifySize(ve1,length(time{j}));
    Pe_srb = ModifySize(Pe_srb,length(time{j}));
    Tc1 = ModifySize(Tc1,length(time{j}));
    Pc_srb = ModifySize(pc_graph,length(time{j})); %psi
    cstar1_srb = ModifySize(cstar1,length(time{j}));
    mdotv_srb{j} = mdot_srb{j}.*ve1_srb;

    while abs(difmp2 - 1) > tol
        tb2(j) = tb2(j) + 1*(1-difmp2);
        time_srb2 = linspace(0,tb2(j),length(pc_graph2)); %sec
        pc_srb2 = ModifySize(pc_graph2,length(time_srb2)); %psi
        cstar2_srb = ModifySize(cstar2,length(time_srb2));
        Pc_srb2 = pc_srb2.*6894.76; %Pa
        mdot_srbB{j} = Pc_srb2.*At_srb2./cstar2_srb;
        intmdot2 = cumtrapz(mdot_srbB{j})*(time_srb2(2)-time_srb2(1));
        mdot_check2 = intmdot2(end);

        difmp2 = mdot_check2/mp2{j};
    end

    time2{j} = linspace(0,time_srb2(end),time_srb2(end)/dt + 1);
    ve2_srb = ModifySize(ve2,length(time2{j}));
    Pe_srb2 = ModifySize(Pe_srb2,length(time2{j}));
    Tc2 = ModifySize(Tc2,length(time2{j}));
    Pc_srb2 = ModifySize(pc_graph2,length(time2{j})); %psi
    cstar2_srb = ModifySize(cstar2,length(time2{j}));
    mdot_srb2{j} = interp1(time_srb2,mdot_srbB{j},time2{j});
    mdotv_srb2{j} = mdot_srb2{j}.*ve2_srb;

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
            F_srb{j}(i+1) = (mdotv_srb{j}(i) + Ae_srb(j).*(Pe_srb(i) - P{j}(i)));
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
            F_srb{j}(i+1) = (mdotv_srb2{j}(k) + Ae_srb2(j).*(Pe_srb2(k) - P{j}(i)));
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