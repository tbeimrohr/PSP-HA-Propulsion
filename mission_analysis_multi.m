clc
clear
%% Mission Analysis for PSP High Altitude Preliminary Study
% 12/20/21
% Thomas Beimrohr
%% Constants
g = 9.81; %m/s2
re = 6378; %km
mu = 3.986*10^5; %km3/s2
z = 100; %km
latitude = 40.425869; %deg
day = 24*60*60; %sec
eta_isp = .84; %could be up to .87
Isp = [200 225 250 275 284.11]; %sec
massfraction = .875;
massfraction2 = .85;
Pc = 6894.757*1000; %pa
cstar = 1765.1; %m/s
At = ((.75/12)^2)*pi*0.092903; %m2
tb = 0;
%% Calculations
%% Ideal Delta V Calculations

dv_100km = sqrt((2*(g/1000)*re*z)/(re+z)); %km/s
surface_velocity = (2*pi*re*cosd(latitude))/day; %km/s
tot_dv_mission_ideal = dv_100km; %km/s
Vc = sqrt(mu/(re+z)) - surface_velocity; %km/s
dv_drag_loss = .4*tot_dv_mission_ideal;
grav_loss = g*30/(1000);
tot_dv_mission = tot_dv_mission_ideal + dv_drag_loss + grav_loss;

scale = 0:.01:1;
dv1 = scale.*tot_dv_mission;
dv2 = (1-scale).*tot_dv_mission;
%% Vehicle Mass

for i = 1:length(Isp)
    MR1(:,i) = exp((dv1.*1000)./(g.*Isp(i)));
    MR2(:,i) = exp((dv2.*1000)./(g.*Isp(i)));

    mpl = 25; %kg
    mp2(:,i) = mpl.*((MR2(:,i)-1)./(MR2(:,i)-((MR2(:,i)-1)./massfraction2))); %kg
    mi2(:,i) = mp2(:,i).*((1-massfraction2)./massfraction2); %kg
    mo2(:,i) = mi2(:,i) + mpl + mp2(:,i); %kg
    mf2(:,i) = mi2(:,i) + mpl; %kg

    mp(:,i) = mo2(:,i).*((MR1(:,i)-1)./(MR1(:,i)-((MR1(:,i)-1)./massfraction))); %kg
    mi(:,i) = mp(:,i).*((1-massfraction)./massfraction); %kg
    mo(:,i) = mi(:,i) + mo2(:,i) + mp(:,i); %kg
    mf(:,i) = mi(:,i) + mo2(:,i); %kg

    tot_mp(:,i) = mp(:,i) + mp2(:,i);
    index = find(tot_mp(:,i) == min(tot_mp(:,i)));
    tot_mp_min(i) = tot_mp(index,i);
    mp_min(i) = mp(index,i);
    mp2_min(i) = mp2(index,i);
    mo_min(i) = mo(index,i);
    mi_min(i) = mi(index,i) + mi2(index,i);
    ideal_spread(i) = scale(index);
end

% Plots
figure(1)
plot(scale,tot_mp)
xlabel('First Stage \DeltaV as a Percentage to Total \DeltaV')
ylabel('Total Rocket Propellant Mass [kg]')
title('Optimal \DeltaV Spread Between Stages')
legend('200s ISP','225s ISP','250s ISP','275s ISP','284s ISP','location','best')
grid on
%hold on
%plot(ideal_spread,tot_mp_min,'rx')


figure(2)
Details = [tot_mp_min', mpl.*ones(length(Isp),1), mi_min', mo_min'];
T = array2table(Details,'VariableNames',{'Total Propellent Mass' 'Payload Mass' 'Total Inert Mass' 'Total Mass'},'RowNames',{'200s ISP';'225s ISP';'250s ISP';'275s ISP';'284s ISP'});
T = table(T,'VariableNames',{'Vechicle Mass in kg with Various ISP Values'}); % Nested table
TString = evalc('disp(T)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
set(gcf,'Position',[100 100 600 200])
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

figure(3)
Details = [tot_dv_mission_ideal,dv_drag_loss,grav_loss,tot_dv_mission];
T = array2table(Details,'VariableNames',{'Ideal','Drag','Gravity','Total'});
T = table(T,'VariableNames',{'\DeltaV in km/s Required for Mission'}); % Nested table
TString = evalc('disp(T)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
set(gcf,'Position',[100 100 400 100])
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
% %% Motor Classification
% Impulse1 = Isp*mp*g;
% Impulse2 = Isp*mp2*g;
% classify = input("Classify the Motor? y = 1 n = 0\n");
% if classify
%     if (Impulse1 > 640.01 & Impulse1 < 1280)
%         classification1 = "J"
%     elseif (Impulse1 > 1208.01 & Impulse1 < 2560)
%         classification1 = "K"
%     elseif (Impulse1 > 2560.01 & Impulse1 < 5120)
%         classification1 = "L"
%     elseif (Impulse1 > 5120.01 & Impulse1 < 10240)
%         classification1 = "M"
%     elseif (Impulse1 > 10240.01 & Impulse1 < 20480)
%         classification1 = "N"
%     elseif (Impulse1 > 20480.01 & Impulse1 < 40960)
%         classification1 = "O"
%     elseif (Impulse1 > 40960.01 & Impulse1 < 81920)
%         classification1 = "P"
%         Impulse1
%     elseif (Impulse1 > 81920.01 & Impulse1 < 163840)
%         classification1 = "Q"
%     elseif (Impulse1 > 163840.01 & Impulse1 < 327680)
%         classification1 = "R"
%     elseif (Impulse1 > 327680.01 & Impulse1 < 655360)
%         classification1 = "s"
%     else
%         classification1 = "BFR"
%     end
%
%     if (Impulse2 > 640.01 & Impulse2 < 1280)
%         classification2 = "J"
%     elseif (Impulse2 > 1208.01 & Impulse2 < 2560)
%         classification2 = "K"
%     elseif (Impulse2 > 2560.01 & Impulse2 < 5120)
%         classification2 = "L"
%     elseif (Impulse2 > 5120.01 & Impulse2 < 10240)
%         classification2 = "M"
%     elseif (Impulse2 > 10240.01 & Impulse2 < 20480)
%         classification2 = "N"
%     elseif (Impulse2 > 20480.01 & Impulse2 < 40960)
%         classification2 = "O"
%         Impulse2
%     elseif (Impulse2 > 40960.01 & Impulse2 < 81920)
%         classification2 = "P"
%     elseif (Impulse2 > 81920.01 & Impulse2 < 163840)
%         classification2 = "Q"
%     elseif (Impulse2 > 163840.01 & Impulse2 < 327680)
%         classification2 = "R"
%     elseif (Impulse2 > 327680.01 & Impulse2 < 655360)
%         classification2 = "s"
%     else
%         classification = "BFR"
%     end
% end
