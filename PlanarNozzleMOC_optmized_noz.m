%%
%  Method of Characteristics Code for 2-D Planar Nozzles
%  ME 510 Spring 2022
%  Aaron Morris
clc
warning('off')
gatherdata = input('Do you want to run characteristic modeling? y=1 n=0\n');
if gatherdata
    clear
    gatherdata = 1;
end

%% User Inputs
%--------------------------------------------------------------------------
% Numerical Properties
%--------------------------------------------------------------------------
%  Number of incident waves
nwave = 25;
%%
if gatherdata
    for runs = 1:length(nwave)
        clc
        fprintf('Run Number %d out of %d...\n',runs,length(nwave))
        %% Define Variables
        %% User Inputs
        
        
        %--------------------------------------------------------------------------
        % Gas Properties
        %--------------------------------------------------------------------------
        %  Gamma
        gamma = 1.4;
        R = 287; %j/kgk
        %  Stagnation pressure [atm] and temperature [K]
        p_stag = 10;
        T_stag = 300;
        rho_stag = convpres(p_stag,'atm','pa')/(R*T_stag);
        
        %--------------------------------------------------------------------------
        % Geometry
        %--------------------------------------------------------------------------
        % Nozzle Shape Selection
        Conical = 0;
        Bell = 1;
        
        %  Nozzle Length
        axial_dist = 1.0;
        %  Throat height
        h_throat = 0.1;
        
        if Conical
            theta_t_deg = 10; %degrees
        elseif Bell
            %  Half-Angle at the throat (in degrees)
            exit_mach = 3;
            [~,nu_opt,~] = flowprandtlmeyer(gamma,exit_mach,'mach');
            theta_t_deg = .5*nu_opt;
        else
            fprintf("Select a Nozzle Shape.\n")
        end
        
        
        
        
        
        %% Calculated Variables
        theta_t = theta_t_deg*pi/180.0;
        
        % Calculates the height at the exit plane for a conical nozzle
        h_exit = 2.0*(tan(theta_t) * axial_dist+0.5*h_throat);
        % Area ratio
        ARatio = h_exit / h_throat;
        
        %% Declared Variables
        % The first index steps through the number of left-running characteristics
        %    that reflect off the axis of symmetry.  It is initially allocated as
        %    "1", but additional "rows" are added as more reflections are computed
        % The second index steps along the nodes of each left-running
        %    characteristic.
        % The number of input waves (nwave) is the number of waves that form at the
        %    expansion fan at the throat, and the +1 is for the wall node
        % The variable n_left_wave is a counter that increases as additional left
        %    running reflections form.
        
        theta = nan(1, nwave(runs)+1); %flow direction at node
        nu = nan(1, nwave(runs)+1); % Prandtl-Meyer function evaluated at node
        x = nan(1, nwave(runs)+1); % x-coordinate of node
        y = nan(1, nwave(runs)+1); % y-coordinate of node
        M = nan(1, nwave(runs)+1); % Mach number
        mu = nan(1, nwave(runs)+1); % Mach Angle
        
        theta_left_char = nan(1, nwave(runs)+1); %flow direction at node
        nu_left_char = nan(1, nwave(runs)+1); % Prandtl-Meyer function evaluated at node
        x_left_char = nan(1, nwave(runs)+1); % x-coordinate of node
        y_left_char = nan(1, nwave(runs)+1); % y-coordinate of node
        M_left_char = nan(1, nwave(runs)+1); % Mach number
        mu_left_char = nan(1, nwave(runs)+1); % Mach Angle
        
        n_left_wave = 1; % This is the number of left running waves
        
        
        %--------------------------------------------------------------------------
        %% Initialize nodes at the throat
        %--------------------------------------------------------------------------
        theta_min_throat = 0.0375*pi/180.0; % Flow direction at the throat (rad),
        %   initally set to a small number.
        
        % delta theta between flow direction
        dtheta = (theta_t - theta_min_throat)/(double(nwave(runs))-1.0);
        
        % The nodes at the throat (corner) are the first row of theta, nu, etc.
        % matrices
        
        % Loop through nodes (starting at #2).  The reason we start at #2 is for
        % book-keeping, so that the number of nodes at the throat equals the number
        % of nodes for each left-running reflection.  After assigning the nodes,
        % copy column 2 to column 1.
        for i = 2: nwave(runs)+1
            theta(1,i) = theta_min_throat + double(i-2)*dtheta;
            nu(1,i) = theta(1,i);
            M(1,i) = InvPrandtlMeyer(gamma, nu(1,i));
            mu(1,i) = asin(1.0/M(1,i));
            x(1,i) = 0.0;
            y(1,i) = h_throat*0.5;
        end
        % Copy one more node to throat (for visualization)
        theta(1,1) = theta(1,2);
        nu(1,1) = nu(1,2);
        M(1,1) = M(1,2);
        mu(1,1) = mu(1,2);
        x(1,1) = x(1,2);
        y(1,1) = y(1,2);
        
        %--------------------------------------------------------------------------
        %% Step through the nodes on the reflected left running wave(s)
        %--------------------------------------------------------------------------
        
        % This is a flag for whether or not to include calculation of the
        % left-running reflection.  This flag turns false when the x-position of
        % the node exceeds the length of the nozzle.
        add_left_reflection = 1;
        
        % Continue adding expansions until the left reflected expansion's
        % x-position exceeds the length of the nozzle
        while add_left_reflection
            
            % Handle the boundary node along the axis of symmetry
            j = n_left_wave;
            % The first reflected wave is a special case where the incident expansion
            % wave is presumed to be linear (as suggested by Anderson).  Therefore the
            % flow angle is a very small (but non zero) value.  This isn't physical,
            % but it is an approximation used by Anderson to maintain the first wave is
            % in fact linear.
            
            % The nodes along the reflected left-running wave a stored in "vectors"
            % for theta_left_char, nu_left_char, etc., and are added to the overall
            % matrix if add_left_reflection is true.
            
            
            % First Step: Calculate the properties at the nodes along the
            % centerline
            
            if j==1 % (Special case for the first reflected wave)
                theta_left_char(1, 1) = theta(j,2);
            else
                % Normal case where the flow direction along the symmetry axis is
                % horizontal
                theta_left_char(1, 1) = 0.0;
            end
            nu_left_char(1, 1) = nu(j,2)+theta(j,2);
            M_left_char(1,1) = InvPrandtlMeyer(gamma,nu_left_char(1,1));
            mu_left_char(1,1) = asin(1.0/M_left_char(1,1));
            mu_avg = 0.5*(mu_left_char(1,1)+mu(j,2));
            theta_avg = 0.5*(theta_left_char(1,1)+theta(j,2));
            x_left_char(1,1) = x(j,2)-y(j,2)/tan(theta_avg-mu_avg);
            y_left_char(1,1) = 0.0;
            
            
            % Second Step: Calculate the internal nodes in non-simple region
            
            for i = 2: nwave(runs)
                % The variables with the "_a" are values for the closest neighbor
                % (upstream) node along the C- characteristic
                theta_a = theta(j,i+1);
                nu_a = nu(j,i+1);
                M_a = M(j,i+1);
                mu_a = mu(j,i+1);
                x_a = x(j,i+1);
                y_a = y(j,i+1);
                
                % The variables with the "_b" are values for the closest neighbor
                % (upstream) node along the C+ characteristic
                theta_b = theta_left_char(1,i-1);
                nu_b = nu_left_char(1,i-1);
                M_b = InvPrandtlMeyer(gamma,nu_b);
                mu_b = asin(1.0/M_b);
                x_b = x_left_char(1,i-1);
                y_b = y_left_char(1,i-1);
                
                % Complete the equations here for theta, nu, M, mu, x, and y
                theta_left_char(1,i) = .5*(theta_a+theta_b) + .5*(nu_a - nu_b);
                nu_left_char(1,i) = .5*(theta_a-theta_b) + .5*(nu_a+nu_b);
                M_left_char(1,i) = InvPrandtlMeyer(gamma,nu_left_char(1,i));
                mu_left_char(1,i) = asin(1.0/M_left_char(1,i));
                
                % Complete the code here to evaluate the position of the nodes
                mac = tan(.5*(theta_a+theta_left_char(1,i)) - .5*(mu_a+mu_left_char(1,i)));
                mbc = tan(.5*(theta_b+theta_left_char(1,i)) + .5*(mu_b+mu_left_char(1,i)));
                x_left_char(1,i) = (mbc*x_b - mac*x_a + y_a - y_b)/(mbc-mac);
                y_left_char(1,i) = y_b + mbc*(x_left_char(1,i) - x_b);
                
                
            end
            % Third Step: Calculate the wall nodes
            
            % -----------------------------------------------------------------
            %  Task 3 - Source code that needs to be modified
            % -----------------------------------------------------------------
            i = nwave(runs)+1;
            % calculate theta, nu, M, and mach angle
            
            % set the angle of flow at the wall
            if Conical
                theta_left_char(1,i) = theta_t;
            elseif Bell
                theta_left_char(1,i) = ((theta_left_char(1,i-1)+theta_b)/2);
            else
                fprintf('Choose a Nozzle Shape.\n')
            end
            % calculate the value of nu at the node
            nu_left_char(1,i) = nu_left_char(1,i-1)+theta_left_char(1,i)-theta_left_char(1,i-1);
            % calculate the Mach number
            M_left_char(1,i) = InvPrandtlMeyer(gamma,nu_left_char(1,i));
            % calculate the Mach angle
            mu_left_char(1,i) = asin(1.0/M_left_char(1,i));
            
            
            % Calculate x and y coordinates of wall node
            % The _a values are the values at the wall node just upstream
            slope_a = tan(0.5*(theta_left_char(1,i)+theta(j,i)));
            x_a = x(j,i);
            y_a = y(j,i);
            % the _b values are the nearest node on the left-running wave
            % approaching the wall
            slope_b = tan(0.5*(theta_left_char(1,i)+theta_left_char(1,i-1))...
                + 0.5*(mu_left_char(1,i)+mu_left_char(1,i-1)));
            x_b = x_left_char(1,i-1);
            y_b = y_left_char(1,i-1);
            
            % calculate x and y position of the node at the wall.
            x_left_char(1,i) = (y_a-y_b+slope_b*x_b-slope_a*x_a)/(slope_b-slope_a);
            y_left_char(1,i) = y_a+slope_a*(x_left_char(1,i)-x_a);
            
            % -------------------------------------------------------------------
            %  End of section of code that should be modified for task 3
            % -------------------------------------------------------------------
            
            % Add the left characteristic to array
            theta = [theta; theta_left_char];
            nu = [nu; nu_left_char];
            M = [M; M_left_char];
            mu = [mu; mu_left_char];
            x = [x; x_left_char];
            y = [y; y_left_char];
            n_left_wave = n_left_wave+1;
            if x_left_char(1,1)>axial_dist
                add_left_reflection = 0;
            end
        end
        
        %--------------------------------------------------------------------------
        % Final Steps: Post-processing and outputing the data
        %--------------------------------------------------------------------------
        Pressure{runs} = p_stag./((1+.5.*(gamma-1).*M.^2).^(gamma/(gamma-1)));
        Temperature{runs} = T_stag./(1+.5.*(gamma-1).*M.^2);
        Son{runs} = sqrt(gamma.*R.*Temperature{runs});
        Vel{runs} = M.*Son{runs};
        Theta{runs} = theta;
        Nu{runs} = nu;
        Mach{runs} = M;
        Mu{runs} = mu;
        X{runs} = x;
        Y{runs} = y;
        stop(runs) = n_left_wave;
    end
    fprintf('Gathering 1D Data...\n')
    %% Post Processing 1D Flow Theory
    
    for i = 1:length(nwave)
        %     eps{i} = Y{i}(:,end)./h_throat;
        wallheight{i} = 2.*tand(theta_t_deg)*X{i} + h_throat;
        eps{i} = (wallheight{i})./h_throat;
        for j = 1:length(eps{i}(1,:))
            for k = 1:length(eps{i}(:,1))
                [Mach_1D{i}(k,j),Temperature_1D_ratio{i}(k,j),Pressure_1D_ratio{i}(k,j),Rho_1D_ratio{i}(k,j),Area_1D_ratio{i}(k,j)] = flowisentropic(gamma,eps{i}(k,j),'sup');
                Temperature_1D{i}(k,j) = Temperature_1D_ratio{i}(k,j).*T_stag;
                Pressure_1D{i}(k,j) = Pressure_1D_ratio{i}(k,j).*p_stag;
                Rho_1D{i}(k,j) = Rho_1D_ratio{i}(k,j).*rho_stag;
                Mach_1D{i}(k,j) = Mach_1D{i}(k,j);
            end
        end
    end
end
fprintf('Plotting Data...\n')
%% Plots
% Plot data along centerline
figure(1)
grid on
Npoint=100; % Number of points to output along the line
xq = linspace(0.0,axial_dist,Npoint);
yq = zeros(1, Npoint);
for i = 1:length(nwave)
    MachQ = griddata(X{i},Y{i},Mach{i},xq,yq); %Interpolate field data to line
    plot(xq,MachQ,'LineWidth',1)
    xlabel('Axial Distance')
    ylabel('Mach Number')
    t2(i) = string(sprintf('%d Waves',nwave(i)));
    leg2 = legend(t2,'location','best');
    title('Mach along the Centerline of the Nozzle')
    hold on
end

%Plot data at exit plane
figure(2)
grid on
Npoint = 50; % Number of points to output along line
xq = axial_dist + zeros(1,Npoint);
yq = linspace(0.0,h_exit, Npoint);
for i = 1:length(nwave)
    MachQ = griddata(X{i},Y{i},Mach{i},xq,yq); %Interpolate field data to line
    plot(MachQ,yq,'LineWidth',1)
    xlabel('Mach Number')
    ylabel('Radial Distance')
    t3(i) = string(sprintf('%d Waves',nwave(i)));
    leg3 = legend(t3,'location','best');
    title('Mach Number along the Exit Plane of the Nozzle')
    hold on
end

figure(3)
grid on
Npoint = 100; % Number of points to output along the line
xq = linspace(0.0,axial_dist,Npoint);
yq = zeros(1, Npoint);
for i = 1:length(nwave)
    PresQ = griddata(X{i},Y{i},Pressure{i},xq,yq); %Interpolate field data to line
    plot(xq,PresQ,'LineWidth',1)
    xlabel('Axial Distance')
    ylabel('Pressure [atm]')
    title('Pressure along the Centerline of the Nozzle')
    t4(i) = string(sprintf('%d Waves',nwave(i)));
    leg4 = legend(t4,'location','best');
    hold on
end

% Plot data at exit plane
figure(4)
grid on
Npoint = 50; % Number of points to output along line
xq = axial_dist + zeros(1,Npoint);
yq = linspace(0.0,h_exit, Npoint);
for i = 1:length(nwave)
    Q = griddata(X{i},Y{i},Pressure{i},xq,yq); %Interpolate field data to line
    plot(Q,yq,'LineWidth',1)
    xlabel('Pressure [atm]')
    ylabel('Radial Distance')
    t3(i) = string(sprintf('%d Waves',nwave(i)));
    leg3 = legend(t3,'location','best');
    title('Pressure along the Exit Plane of the Nozzle')
    hold on
end

figure(5)
subplot(2,1,1)
t = sprintf('Mach with Method of Characteristics with %d Waves',nwave(end));
title(t)
xlabel('Axial Distance')
ylabel('Radial Distance')
hold on
grid on
contourf (x,y,M,500,'LineStyle','none')

%Overlay Nozzle Contour
for j=2:n_left_wave
    plot([x(j-1,nwave(runs)+1);x(j,nwave(runs)+1)],[y(j-1,nwave(runs)+1);y(j,nwave(runs)+1)],...
        'k','LineWidth',3)
end


%Overlay right-running waves

for j=2:n_left_wave
    for i = 1:nwave(runs)
        plot([x(j-1,i+1);x(j,i)],[y(j-1,i+1);y(j,i)],'k')
    end
end
%Overlay left-running waves
for j=2:n_left_wave
    for i = 1: nwave(runs)
        plot([x(j,i);x(j,i+1)],[y(j,i);y(j,i+1)],'k')
    end
end
%Set some limits on the axes
ylim([0, h_exit*.55]);
xlim([0, axial_dist]);
colorbar
colormap('jet')
caxis([0 3.25])
daspect([1,1,1]);
hold off

subplot(2,1,2)
t = sprintf('Mach with Quasi 1D with %d Waves',nwave(end));
title(t)
xlabel('Axial Distance')
ylabel('Radial Distance')
hold on
grid on
contourf (X{end},Y{end},Mach_1D{end},500,'LineStyle','none')

% Overlay Nozzle Contour
for j=2:n_left_wave
    plot([X{end}(j-1,nwave(runs)+1);X{end}(j,nwave(runs)+1)],[Y{end}(j-1,nwave(runs)+1);Y{end}(j,nwave(runs)+1)],...
        'k','LineWidth',3)
end


% Overlay right-running waves

for j=2:n_left_wave
    for i = 1:nwave(runs)
        plot([X{end}(j-1,i+1);X{end}(j,i)],[Y{end}(j-1,i+1);Y{end}(j,i)],'k')
    end
end
% Overlay left-running waves
for j=2:n_left_wave
    for i = 1: nwave(runs)
        plot([X{end}(j,i);X{end}(j,i+1)],[Y{end}(j,i);Y{end}(j,i+1)],'k')
    end
end
% Set some limits on the axes
ylim([0, h_exit*.55]);
xlim([0, axial_dist]);
colorbar
colormap('jet')
caxis([0 3.25])
daspect([1,1,1]);
hold off

figure(6)
grid on
Npoint = 100; % Number of points to output along the line
xq = linspace(0.0,axial_dist,Npoint);
yq = zeros(1, Npoint);
Q = griddata(X{end},Y{end},Mach{end},xq,yq); %Interpolate field data to line
plot(xq,Q,'LineWidth',1)
hold on
Q2 = griddata(X{end},Y{end},Mach_1D{end},xq,yq); %Interpolate field data to line
plot(xq,Q2,'LineWidth',1)
legend('Method of Characteristics','Quasi 1D','location','best')
xlabel('Axial Distance')
ylabel('Mach Number')
title('Comparision of Centerline Mach Between 1D and 2D Methods')

figure(7)
grid on
Npoint = 100; % Number of points to output along the line
xq = linspace(0.0,axial_dist,Npoint);
yq = zeros(1, Npoint);
Q = griddata(X{end},Y{end},Pressure{end},xq,yq); %Interpolate field data to line
plot(xq,Q,'LineWidth',1)
hold on
Q2 = griddata(X{end},Y{end},Pressure_1D{end},xq,yq); %Interpolate field data to line
plot(xq,Q2,'LineWidth',1)
legend('Method of Characteristics','Quasi 1D','location','best')
xlabel('Axial Distance')
ylabel('Pressure [atm]')
title('Comparision of Centerline Pressure Between 1D and 2D Methods')

figure(8)
grid on
Npoint = 50; % Number of points to output along line
xq = axial_dist + zeros(1,Npoint);
yq = linspace(0.0,h_exit, Npoint);
Q = griddata(X{end},Y{end},Mach{end},xq,yq); %Interpolate field data to line
plot(Q,yq,'LineWidth',1)
hold on
Q2 = griddata(X{end},Y{end},Mach_1D{end},xq,yq); %Interpolate field data to line
plot(Q2,yq,'LineWidth',1)
legend('Method of Characteristics','Quasi 1D','location','best')
ylabel('Radial Distance')
xlabel('Mach Number')
title('Comparision of Exit Plane Mach Between 1D and 2D Methods')

figure(9)
grid on
Npoint = 50; % Number of points to output along line
xq = axial_dist + zeros(1,Npoint);
yq = linspace(0.0,h_exit, Npoint);
Q = griddata(X{end},Y{end},Pressure{end},xq,yq); %Interpolate field data to line
plot(Q,yq,'LineWidth',1)
hold on
Q2 = griddata(X{end},Y{end},Pressure_1D{end},xq,yq); %Interpolate field data to line
plot(Q2,yq,'LineWidth',1)
legend('Method of Characteristics','Quasi 1D','location','best')
ylabel('Radial Distance')
xlabel('Pressure [atm]')
title('Comparision of Exit Plane Pressure Between 1D and 2D Methods')

figure(10)
subplot(2,1,1)
t = sprintf('Pressure with the Method of Characteristics with %d Waves',nwave(end));
title(t)
xlabel('Axial Distance')
ylabel('Radial Distance')
hold on
grid on
contourf (X{end},Y{end},Pressure{end},500,'LineStyle','none')

% Overlay Nozzle Contour
for j=2:n_left_wave
    plot([x(j-1,nwave(runs)+1);x(j,nwave(runs)+1)],[y(j-1,nwave(runs)+1);y(j,nwave(runs)+1)],...
        'k','LineWidth',3)
end


% Overlay right-running waves

for j=2:n_left_wave
    for i = 1:nwave(runs)
        plot([x(j-1,i+1);x(j,i)],[y(j-1,i+1);y(j,i)],'k')
    end
end
% Overlay left-running waves
for j=2:n_left_wave
    for i = 1: nwave(runs)
        plot([x(j,i);x(j,i+1)],[y(j,i);y(j,i+1)],'k')
    end
end
% Set some limits on the axes
ylim([0, h_exit*.55]);
xlim([0, axial_dist]);
colorbar
colormap('jet')
caxis([0 5])
daspect([1,1,1]);
hold off
% 
subplot(2,1,2)
t = sprintf('Pressure with Quasi 1D with %d Waves',nwave(end));
title(t)
xlabel('Axial Distance')
ylabel('Radial Distance')
hold on
grid on
contourf (X{end},Y{end},Pressure_1D{end},500,'LineStyle','none')

% Overlay Nozzle Contour
for j=2:n_left_wave
    plot([X{end}(j-1,nwave(runs)+1);X{end}(j,nwave(runs)+1)],[Y{end}(j-1,nwave(runs)+1);Y{end}(j,nwave(runs)+1)],...
        'k','LineWidth',3)
end


% Overlay right-running waves

for j=2:n_left_wave
    for i = 1:nwave(runs)
        plot([X{end}(j-1,i+1);X{end}(j,i)],[Y{end}(j-1,i+1);Y{end}(j,i)],'k')
    end
end
% Overlay left-running waves
for j=2:n_left_wave
    for i = 1: nwave(runs)
        plot([X{end}(j,i);X{end}(j,i+1)],[Y{end}(j,i);Y{end}(j,i+1)],'k')
    end
end
% Set some limits on the axes
ylim([0, h_exit*.55]);
xlim([0, axial_dist]);
colorbar
colormap('jet')
caxis([0 5])
daspect([1,1,1]);
hold off
