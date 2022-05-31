clc
clear
tic

%% Main Body of Code
% TDB 4/4/22
%% Constants

go = 1;
g = 32.174; %ft/s2

% Chamber Geometry Constants
dimension_y = [14.4 72 31.5]; %in refer to homework document "project 1" for labeling vector follows [a b c ...]
dimension_x = [33 365 318.5 384 1390.5]; %in refer to homework document "project 1" for labeling vector follows [A B C ...]
elements = 100;

seg1_axiala = linspace(0,dimension_x(1),elements);
dxseg1a = seg1_axiala(end) - seg1_axiala(end-1);
seg1_axialb = linspace(dimension_x(1)+dxseg1a,dimension_x(2),elements);
seg1_axial = [seg1_axiala seg1_axialb];
gap = linspace(0,1.5,10);
seg2_axial = linspace(0,dimension_x(3),elements)+seg1_axialb(end)+gap(end);
seg3_axial = linspace(0,dimension_x(3),elements)+seg2_axial(end)+gap(end);

seg1web_max = sind(45)*dimension_y(2) - dimension_y(1);
seg1web = linspace(0,seg1web_max,200);
seg2web_max = dimension_y(2) - dimension_y(3);
seg2web = linspace(0,seg2web_max,200);
seg3web = seg2web;
seg4web = seg2web;

% Grain Geometry Constants
inhibited1 = 1;
inhibited2 = 0;
inhibited3 = 0;
inhibited4 = 0;

% Grain Chemistry Constants
CEA_SAVE_FILE = 'CEA_data_new.mat';
p_rho = .063; %llbm/in3
a = .038; %in/(s psi^n)
n = .35;
pc_head_guess = 800; %psi
conv_tol = .01;
pc_range = linspace(75,3000,6*100);
dif_pc = 1;
tol_ss = .01;
% Ru = 1545.349; % (ft*lbf)/(lbmol R)
Ru = 8.31446261815324*1000; %j/mol K

% Nozzle Constants
re = (152.6/2); %in
rt = (53.86/2); %in
Ae = pi*re^2; %in2
At = pi*rt^2; %in2
ep = Ae/At;

% Pc_loss in slots (estimated from AIAA 93-2309)
x_points = [0.15050 0.15860 0.1673 0.18470 0.1985 0.2165 0.2288 0.2549 0.2696 0.28490 0.305 0.3194 0.3389 0.3548 0.3824 0.4091 0.4361 0.4493]';
y_points = ([0.3558823529411761 0.38676470588235246 0.4132352941176465 0.5058823529411764 0.5720588235294115 0.677941176470588 0.7573529411764701 0.9602941176470585 1.0705882352941174 1.216176470588235 1.388235294117647 1.5294117647058822 1.7588235294117647 1.9220588235294116 2.275 2.6455882352941176 3.073529411764706 3.2852941176470587]');

fit_line = fit(x_points,y_points,'exp2');
fit_line_coeff = coeffvalues(fit_line);
x_plot_bounds = linspace(.15,.45,300);
Pt_loss = -.01.*(fit_line_coeff(1)*exp(fit_line_coeff(2).*x_plot_bounds) + fit_line_coeff(3)*exp(fit_line_coeff(4).*x_plot_bounds));

runProgram = input("Do you want to run the balistic model? y=1 n=0\n");
for foldup = 1
if runProgram
    %% Calculations
    if go
        %% Chamber Geometry
        h(1) = ( ( dimension_x(1)^2 + dimension_y(2)^2 - dimension_y(1)^2 )/(2*dimension_x(1)) ); %circular shift in center for the first rounded edge
        R_dim = sqrt( h(1)^2 + dimension_y(1)^2); %Radius for dome
        h(2) = dimension_x(5) - sqrt(R_dim^2 - dimension_y(3)^2);
        seg4_axiala_Stop = h(2) + sqrt(R_dim^2 - dimension_y(2)^2);
        
        seg4_axiala = linspace(seg3_axial(end)+gap(end),seg4_axiala_Stop-gap(end),elements);
        dxseg4a = seg4_axiala(end) - seg4_axiala(end-1);
        seg4_axialb = linspace(seg4_axiala_Stop+gap(end),dimension_x(5),elements+1);
        seg4_axial = [seg4_axiala seg4_axialb(2:end)];
        
        seg1_radiala = sqrt( R_dim^2 - (seg1_axiala - h(1)).^2 );
        seg1_radialb = ones(1,length(seg1_axialb)).*dimension_y(2);
        seg1_radial = [seg1_radiala seg1_radialb];
        seg2_radial = ones(1,length(seg2_axial)).*dimension_y(2);
        seg3_radial = ones(1,length(seg3_axial)).*dimension_y(2);
        seg4_radiala = ones(1,length(seg4_axiala)).*dimension_y(2);
        seg4_radialb = sqrt( R_dim^2 - (seg4_axialb(2:end) -gap(end) - h(2)).^2 );
        seg4_radial = [seg4_radiala seg4_radialb];
        
        radial = [seg1_radiala seg1_radialb seg2_radial seg3_radial seg4_radiala seg4_radialb];
        axial = [seg1_axiala seg1_axialb seg2_axial seg3_axial seg4_axiala seg4_axialb];
        
        dx(1) = 0;
        dy(1) = 0;
        ds(1) = 0;
        for i = 2:length(radial)
            dx(i) = axial(i)-axial(i-1);
            dy(i) = radial(i)-radial(i-1);
            ds(i) = sqrt(dx(i)^2 + dy(i)^2);
        end
        
        dx1 = dx(1:length(seg1_radial));
        dx2 = [0 dx(length(seg1_radial)+2:length(seg2_radial) + length(seg1_radial))];
        dx3 = [0 dx(length(seg2_radial) + length(seg1_radial) + 2:length(seg3_radial) + length(seg1_radial) + length(seg2_radial))];
        dx4 = [0 dx(length(seg3_radial) + length(seg1_radial) + length(seg2_radial) + 2:length(seg4_radial) + length(seg3_radial) + length(seg1_radial) + length(seg2_radial))];
        axial = [cumsum(dx1) (sum(dx1)+gap(end)) cumsum(dx2)+sum(dx1)+gap(end) (sum(dx2)+gap(end)+sum(dx1))  cumsum(dx3)+sum(dx2)+sum(dx1)+gap(end)*2 (sum(dx3)+gap(end)*3+sum(dx2)+sum(dx1))  cumsum(dx4)+sum(dx2)+sum(dx3)+sum(dx1)+gap(end)*3];        %% Grain Geometery
        % Segment 1: Slots
        num_slots = 4;
        for i = 1:length(seg1_radial)
            for j = 1:length(seg1web)
                if j > 1
                    Area1(i,j) = Area1(i,j-1);
                end
                if (dimension_y(1) + seg1web(j)) <= seg1_radial(i)
                    perimeter1(i,j) = num_slots*2*(sqrt(seg1_radial(i)^2 - (dimension_y(1) + seg1web(j)).^2) - (dimension_y(1) + seg1web(j)));
                    theta1(i,j) = asind((dimension_y(1) + seg1web(j))/seg1_radial(i));
                    Area1(i,j) = num_slots*(perimeter1(i,j)/(2*num_slots))*(dimension_y(1)*2 + 2*seg1web(j)) + ((seg1_radial(i)^2)/2)*(deg2rad(theta1(i,j)*2) - sin(deg2rad(theta1(i,j)*2)))*4 + ((dimension_y(1) + seg1web(j))*2)^2;
                    if perimeter1(i,j) <= .0001
                        perimeter1(i,j) = 0;
                        Area1(i,j) = pi*seg1_radial(i)^2;
                    end
                    if i >1
                        AreaG1_2(i,j) = (perimeter1(i,j)*dx1(i)*(dx1(end)/dx1(i)))/(length(dx1)/elements);
                    else
                        AreaG1_2(i,j) = perimeter1(i,j)*dx1(i)*((dx1(end)/dx1(2)))/(length(dx1)/elements);
                    end
                end
                if ~inhibited1
                    AreaG1(i,j) = pi*seg1_radial(i)^2 - Area1(i,j); % Current not written to include an uninhibited seg1
                else
                    AreaG1(i,j) = 0;
                end
                %             if i > 1
                %                 V1(i,j) = dx1(i)*Area1(i,j);
                %             else
                %                 V1(i,j) = dx1(i)*Area1(i,j);
                %             end
            end
            web1_max(i) = find(perimeter1(i,:) == 0,1,'first'); % what the index value is on the segment 1 segment when the propellant burns out
            
        end
        for j = 1:length(seg1web)
            dx1_max(j) = find(perimeter1(:,j) == 0,1,'last');
        end
        V1 = (dx1'.*Area1)./(length(dx1)/elements);
        Ab1 = AreaG1 + AreaG1_2;
        Ab1(Ab1<0) = 0;
        
        %     burnface_index1 = find(Ab1(i,:) == 0,1,'first')
        
        % Segment 2: BATES
        stop = 1;
        k4 = 0;
        for i = 1:length(seg2_radial)
            k = [];
            k3 = 1;
            for j = 1:length(seg2web)
                perimeter2(i,j) = 2*pi*(dimension_y(3) + seg2web(j));
                Area2(i,j) = pi*(dimension_y(3) + seg2web(j))^2; %port area
                AreaG2_2(i,j) = perimeter2(i,j)*dx2(i)/(length(dx2)/elements); %inner core surface area
                
                if ~inhibited2
                    AreaG2(i,j) = 0;
                    AreaL2_2(i,j) = 0;
                    if seg2web(j) >= (seg2_axial(i) - sum(dx1) - gap(end))
                        AreaL2_2(i,j) = ((seg2web(j))*perimeter2(i,j)); %area of part being burned inwards
                        k = find(AreaL2_2(i,:),1,'first');
                        k2 = i;
                    end
                    if i > 1
                        k3 = find(AreaL2_2(i-1,:),1,'first');
                    end
                    AreaG2(k2,k3:(k-1)) = (pi*seg2_radial(k2)^2 - Area2(k2,k3:(k-1))); %area of outer face, the uninhibited side
                    if i < length(seg2_radial) && i > 1 && seg2_axial(i+1) - sum(dx1) - gap(end) >= seg2web_max & stop & ~isempty(k)
                        AreaG2(i,j) = (pi*seg2_radial(i)^2 - Area2(i,j)); %area of outer face, the uninhibited side
                        k4 = 1;
                    end
                else
                    AreaG2(i,j) = 0;
                    AreaL2_2(i,j) = 0;
                end
                %             V2(i,j) = dx2(i)*Area2(i,j);
            end
            if k4
                stop = 0;
            end
            web2_max(i) = length(seg2web);
            dx2_max(i) = 1;
        end
        
        
        AreaG2_3 = AreaG2_2 - AreaL2_2 - flip(AreaL2_2);
        Ab2 = AreaG2 + AreaG2_3 + flip(AreaG2);
        Ab2(Ab2<0) = 0;
        for i = 1:length(seg2_radial)
            Area2(Ab2 == 0) = pi*seg2_radial(i)^2;
        end
        V2 = (dx2'.*Area2)./(length(dx2)/elements);
        
        
        
        % Segment 3: BATES (identical to segment 2)
        perimeter3 = perimeter2;
        Area3 = Area2;
        AreaG3 = AreaG2;
        AreaG3_2 = AreaG2_2;
        Ab3 = Ab2;
        V3 = V2;
        web3_max = web2_max;
        dx3_max = dx2_max;
        
        % Segment 4: BATES
        stop = 1;
        k4 = 0;
        for i = 1:length(seg4_radial)
            k = [];
            k3 = 1;
            for j = 1:length(seg4web)
                if (dimension_y(3) + seg4web(j)) <= round(seg4_radial(i),9)
                    perimeter4(i,j) = 2*pi*(dimension_y(3) + seg4web(j));
                    Area4(i,j) = pi*(dimension_y(3) + seg4web(j))^2;
                    if perimeter4(i,j) <= 0
                        perimeter4(i,j) = 0;
                        Area4(i,j) = pi*seg4_radial(i)^2;
                    end
                    
                    if ~inhibited4
                        %                 AreaG4(i,j) = pi*seg4_radial(i)^2 - Area4(i,j);
                        AreaG4(i,j) = 0;
                        AreaL4_2(i,j) = 0;
                        if seg4web(j) >= (seg4_axial(i) - sum(dx1) - gap(end) - sum(dx2) -gap(end) - sum(dx3) - gap(end))
                            AreaL4_2(i,j) = ((seg4web(j))*perimeter4(i,j)); %area of port being burned inwards
                            k = find(AreaL4_2(i,:),1,'first');
                            k2 = i;
                        end
                        if i > 1
                            k3 = find(AreaL4_2(i-1,:),1,'first');
                        end
                        AreaG4(k2,k3:(k-1)) = (pi*seg4_radial(k2)^2 - Area4(k2,k3:(k-1))); %area of outer face, the uninhibited side
                        if i < length(seg4_radial) && i > 1 && seg4_axial(i+1) - sum(dx1) - gap(end) - sum(dx2) -gap(end) - sum(dx3) - gap(end) >= seg2web_max & stop & ~isempty(k)
                            AreaG4(i,j) = (pi*seg4_radial(i)^2 - Area4(i,j)); %area of outer face, the uninhibited side
                            k4 = 1;
                        end
                    else
                        AreaG4(i,j) = 0;
                    end
                    if i > 1
                        AreaG4_2(i,j) = perimeter4(i,j)*dx4(i)*(dx4(2)/dx4(i))/(length(dx4)/elements);
                    else
                        AreaG4_2(i,j) = perimeter4(i,j)*dx4(i)*(dx4(2)/dx4(2))/(length(dx4)/elements);
                    end
                end
                if Area4(i,j) == 0 & j > 1
                    Area4(i,j) = Area4(i,j-1);
                end
                %             if i > 1
                %                 V4(i,j) = dx4(i)*Area4(i,j);
                %             else
                %                 V4(i,j) = dx4(i)*Area4(i,j);
                %             end
            end
            if k4
                stop = 0;
            end
            if i <= length(dx4)/2 + 1
                web4_max(i) = length(seg4web);
            else
                web4_max(i) = find(perimeter4(i,:) == 0,1,'first');
            end
        end
        for j = 1:length(seg4web)
            if j <= length(seg4_radiala)
                dx4_max(1) = length(seg4_radiala);
            else
                dx4_max(j) = find(perimeter4(:,j) == 0,1,'first')-length(seg4_radiala);
            end
        end
        %     dx4_max(dx4_max == 0) = dx4_max(length(seg4_radiala)+1);
        %     dx4_max = flip(dx4_max);
        
        AreaG4_3 = AreaG4_2 - AreaL4_2;
        Ab4 = AreaG4 + AreaG4_3;
        Ab4(Ab4<0) = 0;
        for i = 1:length(seg4_radial)
            for j = 1:length(seg4web)
                if Ab4(i,j) == 0
                    Area4(i,j) = pi*seg4_radial(i)^2;
                end
            end
        end
        V4 = (dx4'.*Area4)./(length(dx4)/elements);
        
        % Gap Segment
        Area5 = pi*(dimension_y(2))^2;
        Ab5 = 0;
        V5 = gap(end)*Area5.*ones(1,length(seg3web));
        
        
        %% Grain Chemistry
        
        go_cea = input('Do you wish to run CEA? y=1 n=0\n');
        if go_cea
            areyousure = input('Are you sure???\n');
            if areyousure
                data = Grain_Chemistry(pc_range);
                save(CEA_SAVE_FILE, 'data');
            else
                load(CEA_SAVE_FILE)
            end
        else
            load(CEA_SAVE_FILE)
        end
        
        molmass = squeeze(data('mw'));
        molmass = molmass(:,1);
        R = (Ru./molmass); %j/kg-k
        dlv_dlp = squeeze(data('(dlv/dlp)t'));
        dlv_dlp = dlv_dlp(:,1);
        gammas = squeeze(data('gammas'));
        gammas = gammas(:,1);
        gamma = -1.*gammas.*dlv_dlp;
        Tc = squeeze(data('t')); %K
        Tc = Tc(:,1);
        cstar = squeeze(data('cstar')).*3.281; %ft/s
        cstar = cstar(:,1);
        rho = squeeze(data('rho'))./27680; %lbm/in3
        rho = rho(:,1);
        t(1) = 0;
        dt = 1;
        o = 1;
        scale_time = [1 0 0];
        j = 1;
        
        web_check = 1;
        V = [V1' V5' V2' V5' V3' V5' V4']';
        axial_length = sum(dx1) + sum(dx2) + sum(dx3) + sum(dx4) + gap(end)*3;
        tol_dmdt = .01;
        dif_dmdt(1) = 1;
        sign_check(1) = 1;
        flip = 1;
        scale_v = 25;
        constant_scale = scale_v;
        count = 0;
        count_neg = 0;
        count_flip = 1;
        pc_save = pc_head_guess;
        good = 0;
        
        
        while web_check
            Po(1,j) = pc_head_guess(j);
            M(1,j) = 0;
            P(1,j) = Po(1,j);
            mdot(1,j) = 0;
            Tc_seg(1,j) = interp1(pc_range,Tc,Po(1,j));
            gamma_seg(1,j) = interp1(pc_range,gamma,Po(1,j));
            R_seg(1,j) = interp1(pc_range,R,Po(1,j));
            rho_seg(1,j) = interp1(pc_range,rho,Po(1,j));
            cstar_seg(1,j) = interp1(pc_range,cstar,Po(1,j));
            r_seg(1,j) = 0;
            web_seg(1,j) = 0;
            mdot_seg(1,j) = (a*Po(1,j)^n)*p_rho*interp2(seg1web,seg1_axial,Ab1,web_seg(1,j),seg1_axial(1));
            V_seg(1,j) = interp2(seg1web,seg1_axial,V1,web_seg(1,j),seg1_axial(1));
            web_limit(1,j) = 1;
            Aport(1,j) = interp2(seg1web,seg1_axial,Area1,web_seg(1,j),seg1_axial(1));
            web_seg_bit(1,j) = 0;
            Ab_seg(1,j) = interp2(seg1web,seg1_axial,Ab1,web_seg(1,j),seg1_axial(1));
            if j > 1
                dt = sum(rho_seg(:,j-1).*V_seg(:,j-1))/mdot(end,j-1);
            end
            
            % Segment 1
            for i = 2:length(seg1_axial)
                Po(i,j) = P(i-1,j)*( (1 + ((gamma_seg(i-1,j) - 1)/2)*M(i-1,j)^2)^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                r_seg(i,j) = a*Po(i,j)^n;
                web_seg_bit(i,j) = r_seg(i,j)*dt;
                web_seg(i,j) = sum(web_seg_bit(i,:));
                V_seg(i,j) = interp2(seg1web,seg1_axial,V1,web_seg(i,j),seg1_axial(i));
                if isnan(V_seg(i,j))
                    V_seg(i,j) = V_seg(i,j-1);
                end
                Ab_seg(i,j) = interp2(seg1web,seg1_axial,Ab1,web_seg(i,j),seg1_axial(i));
                if isnan(Ab_seg(i,j))
                    Ab_seg(i,j) = 0;
                end
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
                
                if mdot(i,j) == 0
                    Po(i,j) = Po(i-1,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i-1,j)*(M(i-1,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i-1,j)-1)/2)*M(i-1,j)^2 )^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                Aport(i,j) = interp2(seg1web,seg1_axial,Area1,web_seg(i,j),seg1_axial(i));
                if isnan(Aport(i,j))
                    Aport(i,j) = Aport(i,j-1);
                end
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Aport(i,j),R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                
                if mdot(i,j) == 0
                    Po(i,j) = Po(i,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i,j)*(M(i,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i,j)-1)/2)*M(i,j)^2 )^(gamma_seg(i,j)/(gamma_seg(i,j) - 1)));
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
            end
            
            index_seg1 = i;
            
            % Gap 1
            for i = index_seg1+1:index_seg1+1
                Po(i,j) = interp1(x_plot_bounds,Pt_loss,M(i-1,j),"makima",'extrap')*Po(i-1,j) + Po(i-1,j);
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab5)/sum(gap);
                mdot(i,j) = sum(mdot_seg(:,j));
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Area5,R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                Po(i,j) = interp1(x_plot_bounds,Pt_loss,M(i,j),"makima",'extrap')*Po(i-1,j) + Po(i-1,j);
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i,j)-1)/2)*M(i,j)^2 )^(gamma_seg(i,j)/(gamma_seg(i,j) - 1)));
                web_seg(i,j) = 0;
                V_seg(i,j) = V5(1);
            end
            
            index_gap1 = i;
            % Segment 2
            for i = index_gap1+1:length(seg2_axial)+ index_gap1
                
                Po(i,j) = P(i-1,j)*( (1 + ((gamma_seg(i-1,j) - 1)/2)*M(i-1,j)^2)^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                r_seg(i,j) = a*Po(i,j)^n;
                web_seg_bit(i,j) = r_seg(i,j)*dt;
                web_seg(i,j) = sum(web_seg_bit(i,:));
                %
                %             if abs(((web_seg(i,j)-seg2web(web2_max(i- index_gap1)))/seg2web(web2_max(i- index_gap1)))) < .25 && scale_time(1)
                %                 dt = dt/2;
                %                 scale_time(1) = 0;
                %                 scale_time(2) = 1;
                %             elseif abs(((web_seg(i,j)-seg2web(web2_max(i- index_gap1)))/seg2web(web2_max(i- index_gap1)))) < .1 && scale_time(2)
                %                 dt = dt/2;
                %                 scale_time(2) = 0;
                %                 scale_time(3) = 1;
                %             elseif abs(((web_seg(i,j)-seg2web(web2_max(i- index_gap1)))/seg2web(web2_max(i- index_gap1)))) < .0375 && scale_time(3)
                %                 dt = dt/4;
                %                 scale_time(3) = 0;
                %             end
                V_seg(i,j) = interp2(seg2web,seg2_axial,V2,web_seg(i,j),seg2_axial(i - index_gap1));
                if isnan(V_seg(i,j))
                    V_seg(i,j) = V_seg(i,j-1);
                end
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                Ab_seg(i,j) = interp2(seg2web,seg2_axial,Ab2,web_seg(i,j),seg2_axial(i - index_gap1));
                if isnan(Ab_seg(i,j))
                    Ab_seg(i,j) = 0;
                end
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
                if mdot(i,j) == 0
                    Po(i,j) = Po(i-1,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i,j)*(M(i-1,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i-1,j)-1)/2)*M(i-1,j)^2 )^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                Aport(i,j) = interp2(seg2web,seg2_axial,Area2,web_seg(i,j),seg2_axial(i - index_gap1));
                if isnan(Aport(i,j))
                    Aport(i,j) = Aport(i,j-1);
                end
                
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Aport(i,j),R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                if mdot(i,j) == 0
                    Po(i,j) = Po(i,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i,j)*(M(i,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i,j)-1)/2)*M(i,j)^2 )^(gamma_seg(i,j)/(gamma_seg(i,j) - 1)));
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
            end
            
            index_seg2 = i;
            
            % Gap 2
            for i = index_seg2+1:index_seg2+1
                Po(i,j) = interp1(x_plot_bounds,Pt_loss,M(i-1,j),"makima",'extrap')*Po(i-1,j) + Po(i-1,j);
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab5)/sum(gap);
                mdot(i,j) = sum(mdot_seg(:,j));
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Area5,R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                Po(i,j) = interp1(x_plot_bounds,Pt_loss,M(i,j),"makima",'extrap')*Po(i-1,j) + Po(i-1,j);
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i,j)-1)/2)*M(i,j)^2 )^(gamma_seg(i,j)/(gamma_seg(i,j) - 1)));
                web_seg(i,j) = 0;
                V_seg(i,j) = V5(1);
            end
            
            index_gap2 = i;
            % Segment 3
            for i = index_gap2+1:length(seg3_axial)+ index_gap2
                Po(i,j) = P(i-1,j)*( (1 + ((gamma_seg(i-1,j) - 1)/2)*M(i-1,j)^2)^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                r_seg(i,j) = a*Po(i,j)^n;
                web_seg_bit(i,j) = r_seg(i,j)*dt;
                web_seg(i,j) = sum(web_seg_bit(i,:));
                burnface_index3 = find(Ab3(i-index_gap2,:) == 0,1,'first');
                
                Ab_seg(i,j) = interp2(seg3web,seg3_axial,Ab3,web_seg(i,j),seg3_axial(i - index_gap2));
                if isnan(Ab_seg(i,j))
                    Ab_seg(i,j) = 0;
                end
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
                V_seg(i,j) = interp2(seg3web,seg3_axial,V3,web_seg(i,j),seg3_axial(i - index_gap2));
                if isnan(V_seg(i,j))
                    V_seg(i,j) = V_seg(i,j-1);
                end
                if mdot(i,j) == 0
                    Po(i,j) = Po(i-1,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i-1,j)*(M(i-1,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i-1,j)-1)/2)*M(i-1,j)^2 )^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                Aport(i,j) = interp2(seg3web,seg3_axial,Area3,web_seg(i,j),seg3_axial(i - index_gap2));
                if isnan(Aport(i,j))
                    Aport(i,j) = Aport(i,j-1);
                end
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Aport(i,j),R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                
                if mdot(i,j) == 0
                    Po(i,j) = Po(i,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i,j)*(M(i,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i,j)-1)/2)*M(i,j)^2 )^(gamma_seg(i,j)/(gamma_seg(i,j) - 1)));
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
            end
            
            
            index_seg3 = i;
            
            % Gap 3
            for i = index_seg3+1:index_seg3+1
                Po(i,j) = interp1(x_plot_bounds,Pt_loss,M(i-1,j),"makima",'extrap')*Po(i-1,j) + Po(i-1,j);
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab5)/sum(gap);
                mdot(i,j) = sum(mdot_seg(:,j));
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Area5,R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                Po(i,j) = interp1(x_plot_bounds,Pt_loss,M(i,j),"makima",'extrap')*Po(i-1,j) + Po(i-1,j);
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i-1,j)-1)/2)*M(i-1,j)^2 )^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                web_seg(i,j) = 0;
                V_seg(i,j) = V5(1);
            end
            
            index_gap3 = i;
            % Segment 4
            for i = index_gap3+1:(length(seg4_axial)+ index_gap3)
                
                Po(i,j) = P(i-1,j)*( (1 + ((gamma_seg(i-1,j) - 1)/2)*M(i-1,j)^2)^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                r_seg(i,j) = a*Po(i,j)^n;
                web_seg_bit(i,j) = r_seg(i,j)*dt;
                web_seg(i,j) = sum(web_seg_bit(i,:));
                %             burnface_index4 = find(Ab4(i-index_gap3,:) == 0,1,'first');
                %             if ~isempty(burnface_index4)
                %                 web_seg_bit(i,j) = web_seg_bit(i,j)+ r_seg(i,j)*dt;
                %                 web_seg(i,j) = sum(web_seg_bit(i,:));
                %             end
                %
                Ab_seg(i,j) = interp2(seg4web,seg4_axial,Ab4,web_seg(i,j),seg4_axial(i - index_gap3));
                if isnan(Ab_seg(i,j))
                    Ab_seg(i,j) = 0;
                end
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
                V_seg(i,j) = interp2(seg4web,seg4_axial,V4,web_seg(i,j),seg4_axial(i - index_gap3));
                if isnan(V_seg(i,j))
                    V_seg(i,j) = V_seg(i,j-1);
                end
                if mdot(i,j) == 0
                    Po(i,j) = Po(i-1,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i-1,j)*(M(i-1,j)^2)*(mdot_seg(i,j)/mdot(i,j)));
                end
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i-1,j)-1)/2)*M(i-1,j)^2 )^(gamma_seg(i-1,j)/(gamma_seg(i-1,j) - 1)));
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                Aport(i,j) = interp2(seg4web,seg4_axial,Area4,web_seg(i,j),seg4_axial(i - index_gap3));
                if isnan(Aport(i,j))
                    Aport(i,j) = Aport(i,j-1);
                end
                M(i,j) = NewtonsMethod_P1(mdot_seg(i,j),mdot(i,j),Po(i,j),Aport(i,j),R_seg(i,j),Tc_seg(i,j),gamma_seg(i,j),g,M(i-1,j));
                
                if mdot(i,j) == 0
                    Po(i,j) = Po(i,j);
                else
                    Po(i,j) = Po(i-1,j)/(1 + gamma_seg(i,j)*(M(i,j)^2)*(mdot_seg(i,j)/mdot(i,j)) );
                end
                
                Tc_seg(i,j) = interp1(pc_range,Tc,Po(i,j));
                gamma_seg(i,j) = interp1(pc_range,gamma,Po(i,j));
                R_seg(i,j) = interp1(pc_range,R,Po(i,j));
                rho_seg(i,j) = interp1(pc_range,rho,Po(i,j));
                cstar_seg(i,j) = interp1(pc_range,cstar,Po(i,j));
                P(i,j) = Po(i,j)/((1 + ((gamma_seg(i,j)-1)/2)*M(i,j)^2 )^(gamma_seg(i,j)/(gamma_seg(i,j) - 1)));
                mdot_seg(i,j) = ((a*Po(i,j)^n)*p_rho*Ab_seg(i,j));
                mdot(i,j) = sum(mdot_seg(:,j));
                
                
            end
            index_seg4 = i;
            mass_in(j) = sum(rho_seg(:,j).*V_seg(:,j));
            
            
            if pc_head_guess(j) < 200
                web_check = 0;
            end
            if j > 1
                count = count + 1;
                % convergence
                dm_dt_acc(j) = sum(rho_seg(:,j).*V_seg(:,j))*(((Po(end,j)-Po(end,j-1))/dt)*(1/Po(end,j)) + ((sum(V_seg(:,j))-sum(V_seg(:,j-1)))/dt)*(1/sum(V_seg(:,j))));
                dm_dt_exp(j) = mdot(end,j) - (g*Po(end,j)*At)/cstar_seg(end,j);
                dif_dmdt(j) = dm_dt_acc(j)/dm_dt_exp(j);
                sign_check(j) = sign(dif_dmdt(j));
                %             dif_dmdt(j)
                if dm_dt_acc(j) > (1+conv_tol)*dm_dt_exp(j)
                    guess_high = pc_head_guess(j);
                    pc_head_guess(j) = (guess_high+guess_low)/2;
                end
                if dm_dt_acc(j) < (1-conv_tol)*dm_dt_exp(j)
                    guess_low = pc_head_guess(j);
                    pc_head_guess(j) = (guess_high+guess_low)/2;
                end
                
                if sign_check(j) < 0
                    dm_dt_acc_temp = abs(dm_dt_acc(j));
                    dm_dt_exp_temp = abs(dm_dt_exp(j));
                    if (1-conv_tol)*dm_dt_exp(j) <= dm_dt_acc_temp && dm_dt_acc_temp <= (1+conv_tol)*dm_dt_exp(j)
                        t(j+1) = t(j) + dt;
                        t(j+1)
                        good = 1;
                        pc_head_guess(j) = pc_head_guess(j-1);
                        guess_low = .7*pc_head_guess(j);
                        guess_high = 1.1*pc_head_guess(j);
                    end
                end
                
                
                %             dif_dmdt(j) = 1;
                if abs(abs(dif_dmdt(j)) - 1) < tol_dmdt || abs(guess_high - guess_low)/2 < 1 || good
                    t(j+1) = t(j) + dt;
                    t(j+1)
                    j = j + 1;
                    good = 0;
                    pc_head_guess(j) = pc_head_guess(j-1);
                    guess_low = .7*pc_head_guess(j);
                    guess_high = 1.1*pc_head_guess(j);
                else
                    web_seg(:,j) = [];
                    Tc_seg(:,j) = [];
                    gamma_seg(:,j) = [];
                    R_seg(:,j) = [];
                    rho_seg(:,j) = [];
                    cstar_seg(:,j) = [];
                    M(:,j) = [];
                    P(:,j) = [];
                    mdot_seg(:,j) = [];
                    mdot(:,j) = [];
                    V_seg(:,j) = [];
                    
                end
            else
                t(j+1) = t(j) + dt;
                j = j + 1;
                pc_head_guess(j) = pc_head_guess(j-1);
                guess_low = .7*pc_head_guess(j);
                guess_high = 1.1*pc_head_guess(j);
            end
        end
    end
else
    load('partb.mat')
end
time = t;
load('trajectory.mat')
Mach_exhaust = squeeze(data('mach'));
Mach_exhaust = Mach_exhaust(:,end);
son_exhaust = squeeze(data('son'));
son_exhaust = son_exhaust(:,end);
Pressure_exhaust = squeeze(data('p'));
Pressure_exhaust = Pressure_exhaust(:,end);

for i = 1:length(time)-1
    index_flight = interp1(t,1:1:length(t),time(i),'nearest');
    [Ta,Sona,Pa(i),Rhoa] = atmoscoesa(alt(index_flight));
    Pressure_exhaust_flight(i) = convpres(interp1(pc_range,Pressure_exhaust,Po(1,i)),'Pa','psi');
    Mach_exhaust_flight(i) = interp1(pc_range,Mach_exhaust,Po(1,i));
    son_exhaust_flight(i) = interp1(pc_range,son_exhaust,Po(1,i));
    ve(i) = Mach_exhaust_flight(i)*son_exhaust_flight(i);
    Pa(i) = convpres(Pa(i),'Pa','psi');
    if time(i) == 1
        Thrust(i) = 0;
    else
        Thrust(i) = 3.281*ve(i)*mdot(end,i)/g + Ae*(Pressure_exhaust_flight(i) - Pa(i));
    end
end

end
%% Plots
figure(1)
[wdata,adata] = meshgrid(seg1web,seg1_axial(1:100));
perimeter1graph1 = interp2(wdata,adata,perimeter1(1:100,:),seg1web,10);
perimeter1graph2 = interp2(wdata,adata,perimeter1(1:100,:),seg1web,20);
perimeter1graph3 = interp2(wdata,adata,perimeter1(1:100,:),seg1web,33);

plot(seg1web(perimeter1graph1>0),perimeter1graph1(perimeter1graph1>0))
hold on
plot(seg1web(perimeter1graph2>0),perimeter1graph2(perimeter1graph2>0))
plot(seg1web(perimeter1graph3>0),perimeter1graph3(perimeter1graph3>0))
plot(seg2web,perimeter2(1,:))
legend('Element at 10in','Element at 20in','Element at 40in','Element at 400in','location','best')
xlabel('Web Distance [in]')
ylabel('Perimeter of Grain [in]')
title('Perimeter of Exposed Grain as a Function of Web')
grid on

figure(2)
[wdata2,adata2] = meshgrid(seg1web,seg1_axial(1:100));
area2graph1 = interp2(wdata2,adata2,Area1(1:100,:),seg1web,10);
area2graph2 = interp2(wdata2,adata2,Area1(1:100,:),seg1web,20);
area2graph3 = interp2(wdata2,adata2,Area1(1:100,:),seg1web,33);

plot(seg1web(area2graph1>0),area2graph1(area2graph1>0))
hold on
plot(seg1web(area2graph2>0),area2graph2(area2graph2>0))
plot(seg1web(area2graph3>0),area2graph3(area2graph3>0))
plot(seg2web,Area2(1,:))
legend('Element at 10in','Element at 20in','Element at 40in','Element at 400in','location','best')
xlabel('Web Distance [in]')
ylabel('Area of Port [in^2]')
title('Area of Port as a Function of Web')
grid on

figure(3)
plot(time,Po(1,:))
hold on
plot(time,Po(end,:))
legend('Head End','Aft End','location','best')
xlabel('Time of Burn [sec]')
ylabel('Stagnation Pressure [psi]')
title('Chamber Pressure History')
grid on

figure(4)
plot(time,[0 mdot(end,:)])
xlabel('Time of Burn [sec]')
ylabel('Exit Mass Flow from Chamber [lbm/s]')
title('Chamber Exit Mass Flow History')
grid on

figure(5)
plot(axial,Po(:,interp1(time,1:1:length(time),0,'nearest')))
hold on
plot(axial,Po(:,interp1(time,1:1:length(time),10,'nearest')))
plot(axial,Po(:,interp1(time,1:1:length(time),60,'nearest')))
plot(axial,Po(:,interp1(time,1:1:length(time),72,'nearest')))
plot(axial,P(:,interp1(time,1:1:length(time),0,'nearest')))
plot(axial,P(:,interp1(time,1:1:length(time),10,'nearest')))
plot(axial,P(:,interp1(time,1:1:length(time),60,'nearest')))
plot(axial,P(:,interp1(time,1:1:length(time),72,'nearest')))
legend('Stagnation @ 0s','Stagnation @ 10s','Stagnation @ 60s','Stagnation @ 72s','Static @ 0s','Static @ 10s','Static @ 60s','Static @ 72s','location','best')
xlabel('Axial Distance [in]')
ylabel('Pressure [psi]')
title('Pressure as a Function of Axial Distance Through Motor')
grid on

figure(6)
plot(axial,M(:,interp1(time,1:1:length(time),0,'nearest')))
hold on
plot(axial,M(:,interp1(time,1:1:length(time),10,'nearest')))
plot(axial,M(:,interp1(time,1:1:length(time),60,'nearest')))
plot(axial,M(:,interp1(time,1:1:length(time),72,'nearest')))
legend('0 seconds','10 seconds','60 seconds','72 seconds','location','best')
xlabel('Axial Distance [in]')
ylabel('Mach')
title('Mach as a Function of Axial Distance Through Motor')
grid on

figure(7)
plot(axial,mdot(:,interp1(time,1:1:length(time),0,'nearest')))
hold on
plot(axial,mdot(:,interp1(time,1:1:length(time),10,'nearest')))
plot(axial,mdot(:,interp1(time,1:1:length(time),60,'nearest')))
plot(axial,mdot(:,interp1(time,1:1:length(time),72,'nearest')))
legend('0 seconds','10 seconds','60 seconds','72 seconds','location','best')
xlabel('Axial Distance [in]')
ylabel('Mass Flow Rate [lbm/s]')
title('Mass Flow as a Function of Axial Distance Through Motor')
grid on

figure(8)
plot(time,[0 Thrust]./(10^6))
xlabel('Time of Burn [sec]')
ylabel('Thrust [10^6 lbf]')
title('Thrust of the Rocket During Flight')
grid on
toc