% Mattia Poinelli
% JPL, december 2018
%
% Tracking propagation rates

%clear all,  clc, close all
% profile on
pizzi
%% LOADING DATA

% Parameters
run Europa_physics.m
run Europa_LEFM_parameters.m

% Europa raster
[raster] = imread('Europa_raster_coarse.jpg');

% Constants
Lambda = (mu/eta)/n;
Z = 0.5 * (n^2*R*mu)/(g*sqrt(1+Lambda^2));

Delta = T_ns/(4*pi*(eta/mu));
C = 0.5 * (n^2*R*mu)/(g*sqrt(1+Delta^2));

COMP = 1/E*[1,  -nu,      0 ;
            -nu,   1,      0 ;
             0,    0, 2+2*nu]; % compliance matrix

%% GRID EVALUATION & STRESS FUNCTIONS
phi      = [0:10:360]; phi_rad = deg2rad(phi);          % longitudes
theta    = 90 - [0:10:180]; theta_rad = deg2rad(theta); % latitudes

% Irwin solution (stress Control)
F = @(a,b)  1 + 0.128.*(a/b) - 0.288.*(a/b).^2 + 1.525.*(a/b).^3;
dF = @(a,b) 0.128/b  - 0.576.*a/b^2 + 4.575.* a.^2/b^3;

% Tada solution  (displacement control)
V = @(a,b)  -0.071-0.535*(a/b)+0.169*(a/b)^2-0.09*(a/b)^3+0.02*(a/b)^4-1.071*(b/a)*log(1-a/b);
dV = @(a,b) -0.535*(1/b)+0.338*(a/b^2)-0.27*(a^2/b^3)+0.08*(a^3/b^4)+1.071*(b/a)*(log(1-a/b)/a+1/(b-a));

% Stress functions
sigma_theta = @(t,cotheta,phi) Z .* (-6.*e.*beta_theta_theta(cotheta,0,h_d,l_d).*cos(n.*t+atan(Lambda)) + ...
    e.*beta_theta_theta(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* sin (n.*t + atan(Lambda)) + ...
    3.*cos(2.*phi).*cos(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_theta_theta(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*sin(w + n.*t + atan(Lambda))));% + ... % beginning of secular parts
    %C .* alpha_theta_theta (cotheta,2,h_s,l_s) .* ...
    %cos(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

sigma_phi = @(t,cotheta,phi) Z .* (-6.*e.*beta_phi_phi(cotheta,0,h_d,l_d).*cos(n.*t+atan(Lambda)) + ...
    e.*beta_phi_phi(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* sin (n.*t + atan(Lambda)) + ...
    3.*cos(2.*phi).*cos(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_phi_phi(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*sin(w + n.*t + atan(Lambda))));% + ...
    %C .* alpha_phi_phi (cotheta,2,h_s,l_s) .* ... % beginning of secular parts
    %cos(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

tau = @(t,cotheta,phi) Z.* (2.*e.*beta_theta_phi(cotheta,2,h_d,l_d) .* ...
    (4.* cos(2.*phi).* sin(n.*t + atan(Lambda)) - ...
    3.* sin(2.*phi).*cos(n.*t+atan(Lambda))) + ...
    4 .* cos(epsilon).*sin(epsilon).* beta_theta_phi(cotheta,1,h_d,l_d).* ...
    sin(phi).*sin(w + n.*t + atan(Lambda)));% - ... % beginning of secular parts
    %C .* alpha_theta_phi (cotheta,2,h_s,l_s) .* ...
    %sin(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

% Temporal derivative of Stress function
d_sigma_theta = @(t,cotheta,phi) n.* Z .* (6.*e.*beta_theta_theta(cotheta,0,h_d,l_d).*sin(n.*t+atan(Lambda)) + ...
    e.*beta_theta_theta(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* cos (n.*t + atan(Lambda)) - ...
    3.*cos(2.*phi).*sin(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_theta_theta(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*cos(w + n.*t + atan(Lambda))));% - ... % beginning of secular parts
    %(4.*pi./T_ns).*C .* alpha_theta_theta (cotheta,2,h_s,l_s) .* ...
    %sin(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

d_sigma_phi = @(t,cotheta,phi) n.*Z .* (6.*e.*beta_phi_phi(cotheta,0,h_d,l_d).*sin(n.*t+atan(Lambda)) + ...
    e.*beta_phi_phi(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* cos (n.*t + atan(Lambda)) - ...
    3.*cos(2.*phi).*sin(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_phi_phi(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*cos(w + n.*t + atan(Lambda))));% + ...
    %(4.*pi./T_ns).*C .* alpha_phi_phi (cotheta,2,h_s,l_s) .* ... % beginning of secular parts
    %sin(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

d_tau = @(t,cotheta,phi) n.*Z.* (2.*e.*beta_theta_phi(cotheta,2,h_d,l_d) .* ...
    (4.* cos(2.*phi).* cos(n.*t + atan(Lambda)) + ...
    3.* sin(2.*phi).*sin(n.*t+atan(Lambda))) + ...
    4 .* cos(epsilon).*sin(epsilon).* beta_theta_phi(cotheta,1,h_d,l_d).* ...
    sin(phi).*cos(w + n.*t + atan(Lambda)));% - ... % beginning of secular parts
    %(4.*pi./T_ns).*C .* alpha_theta_phi (cotheta,2,h_s,l_s) .* ...
    %cos(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

%% MAPPING SETTINGS

% Features loading
S = shaperead('Cycloids_smoothed.shp');

%depth (1,:) =  max(max(depth)); depth (end,: ) =  max(max(depth));
load('Critical_depths_file.mat')
depth = depth(:,:,1:end-2);
Mean_depth = mean(depth,3);

phi      = [0:1:360]; phi_rad = deg2rad(phi);          % longitudes
theta    = 90 - [0:1:180]; theta_rad = deg2rad(theta); % latitudes

%% LEFM LOOP
for q = 1:6  %selection of a the feature from the shape file
    disp('Beginning simulation crack'), q
    time = 0;
    TIME = [];
    
    % avoid NaN values, xx = longitudes and yy = latitudes
    xx = S(q).X; xx = xx(~isnan(xx));
    yy = S(q).Y; yy = yy(~isnan(yy));
    END = length(xx);
    Perc_completion (q) = 100;
    Length_crevasse(q) = 0;
    
    CRACK(q).XX = xx;
    CRACK(q).YY = yy;
    
    % loop on the feature
    for i = 1:END
        
        START = time;
        check = 1;
        Standby = 0;
        %% Rotation matrices
        if i == length(xx)
            % last node of the grid has the same charateristic of the
            % previous node to avoid Matlab unconsistencies
            ROT(:,:,i) =  ROT(:,:,i-1);
            
            Tang_ve(1,i) = Tang_ve(1,i-1);
            Tang_ve(2,i) = Tang_ve(2,i-1);
            Norm_ve(1,i) = Norm_ve(1,i-1);
            Norm_ve(2,i) = Norm_ve(2,i-1);
            
            angular_distance (i) = angular_distance(i-1);
            Normal_stress    (i) = Normal_stress   (i-1);
            d_Normal_stress  (i) = d_Normal_stress (i-1);
            Shear_stress     (i) = Shear_stress    (i-1);
            d_Shear_stress   (i) = d_Shear_stress  (i-1);
            
            % strain rate
            epsilon_dot(i) = epsilon_dot (i-1);
            epsilon(i)     = epsilon     (i-1);
            
            theta_Mohr (i) = theta_Mohr(i-1);
            
        elseif yy(i+1) > yy(i) % 'ascending' feature
            % evaluation of useful angles a,b
            b_ang(i) = deg2rad(xx(i+1) - xx(i));
            a_ang(i) = deg2rad(yy(i+1) - yy(i));
            
            % Side_Angle_Side problem, A is the rotation angle, and the
            % last output c gives the length of the segment
            [A,~,c] = Side_Angle_Side (a_ang(i),b_ang(i));
            
            ROT(:,:,i) = [cos(A),-sin(A);
                sin(A),cos(A)]; % counter clock-wise rotation (positive angle)
            
            % Tangent vector
            Tang_ve(1,i) = ROT(1,1,i);
            Tang_ve(2,i) = ROT(2,1,i);
            
            % Normal vector
            Norm_ve(1,i) = ROT(1,2,i);
            Norm_ve(2,i) = ROT(2,2,i);
            
            angular_distance (i) = c;
            
            theta_Mohr(i) = pi - A;
            CRACK(q).ASCENDING(i) = 1;
            
        elseif yy(i+1) < yy(i)% 'descending' feature
            % evaluation of useful angles a,b
            b_ang(i) = deg2rad(xx(i+1) - xx(i));
            a_ang(i) = abs(deg2rad(yy(i+1) - yy(i)));
            
            % Side_Angle_Side problem, A is the rotation angle, and the
            % last output c gives the length of the segment
            [A,~,c] = Side_Angle_Side (a_ang(i),b_ang(i));
            
            ROT(:,:,i) = [cos(-A),-sin(-A);
                sin(-A),cos(-A)]; % clock-wise rotation (negative angle)
            
            % Tangent vector
            Tang_ve(1,i) = ROT(1,1,i);
            Tang_ve(2,i) = ROT(2,1,i);
            
            % Normal vector
            Norm_ve(1,i) = ROT(1,2,i);
            Norm_ve(2,i) = ROT(2,2,i);
            
            angular_distance (i) = c;
            
            theta_Mohr (i) = A;
            CRACK(q).ASCENDING(i) = 0;
        end
        
        % Crevasse length
        Length_crevasse(q) = Length_crevasse(q) + angular_distance (i)*R;
        
        %% TIME LOOP FOR FINDING ACTIVITY
        while check > 0
            
            a = angular_distance(i)*R;
            
            % Diurnal and secular stress evaluation as Jara-Orue, 2011
            SIGMA_COTHETA = sigma_theta (time,deg2rad(90 - yy(1:i)),deg2rad(xx(1:i)));
            SIGMA_PHI     = sigma_phi   (time,deg2rad(90 - yy(1:i)),deg2rad(xx(1:i)));
            TAU           = tau         (time,deg2rad(90 - yy(1:i)),deg2rad(xx(1:i)));
            
            % Diurnal and secular stress Derivative
            D_SIGMA_COTHETA = d_sigma_theta (time,deg2rad(90 - yy(1:i)),deg2rad(xx(1:i)));
            D_SIGMA_PHI     = d_sigma_phi   (time,deg2rad(90 - yy(1:i)),deg2rad(xx(1:i)));
            D_TAU           = d_tau         (time,deg2rad(90 - yy(1:i)),deg2rad(xx(1:i)));
            
            % Cauchy Stress Tensor for every point of the Feature
            T(:,:,i) = [SIGMA_PHI(end),TAU(end);
                TAU(end),SIGMA_COTHETA(end)];      % [Pa]
            % Derivative Cauchy Stress Tensor for every point of the Feature
            d_T(:,:,i) = [D_SIGMA_PHI(end),D_TAU(end);
                D_TAU(end),D_SIGMA_COTHETA(end)];  % [Pa/s]
            
            % Strain Tensor
            STRAIN_VECT = COMP*[SIGMA_PHI(end),SIGMA_COTHETA(end),TAU(end)]';       % [~]
            STRAIN(:,:,i)= [STRAIN_VECT(1),STRAIN_VECT(3)/2;
                STRAIN_VECT(3)/2,STRAIN_VECT(2)];
            
            % Strain Rate Tensor
            d_STRAIN_VECT = COMP*[D_SIGMA_PHI(end),D_SIGMA_COTHETA(end),D_TAU(end)]'; % [1/s]
            d_STRAIN(:,:,i) = [d_STRAIN_VECT(1),d_STRAIN_VECT(3)/2;
                d_STRAIN_VECT(3)/2,d_STRAIN_VECT(2)];
            
            clear STRAIN_VECT d_STRAIN_VECT
            
            %% Rotation of tensors to local coordinates
            
            % Stress/ strain tensor & derivative  rotated to LOCAL coordinates
            T_R(:,:,i)         = ROT(:,:,i) * T       (:,:,i) * ROT(:,:,i)';
            STRAIN_R (:,:,i)   = ROT(:,:,i) * STRAIN  (:,:,i) * ROT(:,:,i)';
            d_T_R(:,:,i)       = ROT(:,:,i) * d_T     (:,:,i) * ROT(:,:,i)';
            d_STRAIN_R (:,:,i) = ROT(:,:,i) * d_STRAIN(:,:,i) * ROT(:,:,i)';
            
            % Normal and Shear Stress & Derivatives
            Normal_stress   (i) = T_R(1,1,i);
            d_Normal_stress (i) = d_T_R(1,2,i);
            Shear_stress    (i) = T_R(1,2,i);
            d_Shear_stress  (i) = d_T_R(1,2,i);
            
            % Normal displacement rate
            epsilon_dot     (i) = d_STRAIN_R(1,1,i);
            epsilon(i) = STRAIN_R(1,1,i);
            
            % Normal stress for the rest of the feature at this Time
            % Mohr Circles
            Sigma_N_MOHR   (1:i) = 1/2 .* (SIGMA_COTHETA+SIGMA_PHI) + ...
                1/2 .* (SIGMA_PHI - SIGMA_COTHETA) .* cos(2.*theta_Mohr(1:i)) + ...
                TAU .* sin (2.*theta_Mohr(1:i));
      
            if i ~= 1
                for j = i-1:-1:1
                    if Sigma_N_MOHR(j) > 0 && a < 150e3 % limitation at a = 150 km;
                        a = a + angular_distance(j)*R;
                    else
                        break
                    end
                end
                
            end
            % geometrical factor as Tada, 2000
            b = a/0.7;
            
            % Stress Intensity factor
            K_I (i) = Normal_stress(i).*sqrt(pi*a)*F(a,b);
            
            % opening diplacement and rate (from elasticity theory)
            delta (i) = 4* Normal_stress(i)*a*V(a,b)/E;
            opening_rate(i) = delta (i) * epsilon_dot(i); % approximation
            segment_open(i) = a;
            % Propagation rate from the derivative of displacement control
            V_rate(i) =  (E/4*opening_rate(i) - d_Normal_stress(i)*a*V(a,b))/...
                (Normal_stress(i)*V(a,b) + Normal_stress(i)*a*dV(a,b));
            
            % Mohr normal stress (check)
            Normal_stress_M (i) = Sigma_N_MOHR(i);
            
            %% PROPAGATION CHECK
                            
            if K_I(i) > TOU*1000 && V_rate(i) > 0 
                % propagation
                check = -1;
                
                clear a
            elseif time - START >= 2* T_europa
                STRING = strcat('Feature :',num2str(q),' formed at :',num2str(floor(i/END*100)),'%, node',num2str(i));
                disp(STRING)
                Perc_completion (q) = floor(i/END*100);
                
                % Ten cycles have passed, Feature has stopped to propagate
                check = -1;
                END = i;
                
                % No propagation
                delta (i) = 0;
                opening_rate(i) = 0;
                V_rate(i) = 0;

            else
                % not propagation, time passes 
                check = 1;
                time = time + 900;
                clear a 
                Standby = Standby + 900;
                
            end
            
        end
        CRACK(q).STANDBY(i) = Standby;
        % time step for finalising the shaped segment
        time_step(i) = angular_distance(i)*R/2/V_rate(i);
        Segment(i) = angular_distance(i)*R;
        
        time = time + time_step(i);
        TIME = [TIME,time];
       
        if i == END % End of the feature
            disp('End simulation')
            CRACK(q).PROPAGATION_RATES = V_rate;
            CRACK(q).OPENING_RATE = opening_rate;
            CRACK(q).WIDTH = delta;
            CRACK(q).LENGTH = Length_crevasse(q);
            CRACK(q).MEAN_SEGMENT = mean(angular_distance)*R;
            CRACK(q).SEGMENT = Segment;
            CRACK(q).TIME = TIME;
            CRACK(q).EPSILON_DOT = epsilon_dot;
            CRACK(q).SEGMENTOPEN = segment_open;
            CRACK(q).STRESS = Normal_stress;
            CRACK(q).STRESS_M = Normal_stress_M;
            XX = xx;
            YY = yy;
            
%             % Problem with Phoenix linea
%             if q ==24
%                time = TIME(end)-TIME(1);
%             end
            
            %LEGEND{q} = num2str(S(q).Name);
%             Cycles_to_build(q) = time/T_europa;
%             Max_opening_width(q) = max(delta);
%             Max_opening_rate (q) = max(opening_rate); % m/s
%             Max_prop_rate(q) = max(V_rate);
%             Mean_prop_rate (q) = a/time;%*3600/1000; 
                       
            clear a delta opening_rate time_step V_rate epsilon_dot epsilon
            clear TIME Segment segment_open Normal_stress
            clear SIGMA_COTHETA SIGMA_PHI TAU STRAIN
            clear D_SIGMA_COTHETA D_SIGMA_PHI D_TAU d_STRAIN
            clear xx yy
            break
        end
        
    end
    
end

% %% HISTOGRAM PLOT
% ACTIVITY = [CRACK.PROPAGATION_RATES];
% figure('Units','inches','Position',[0 0 3.5 3.5],'PaperPositionMode','auto')
% hold on, grid on
% histogram(log10(ACTIVITY))
% xl1 = xlabel('log(Propagation Rate [m/s])');
% yl1 = ylabel('Number of Nodes');
% 
% set(xl1,'Fontsize',9,'Interpreter','Latex')
% set(yl1,'Fontsize',9,'Interpreter','Latex')
% %set(gca,'xscale','log')

%% PROPAGATION RATE PLOT
%figure('Units','inches','Position',[0 0 7.25 7.25],'PaperPositionMode','auto'), hold on
figure(1),hold on

% C1
subplot(2,2,1),hold on, grid on, k = 1;
xx = S(k).X; xx = xx(~isnan(xx));
yy = S(k).Y; yy = yy(~isnan(yy));

plot(CRACK(k).XX(1:length(CRACK(k).PROPAGATION_RATES)),CRACK(k).PROPAGATION_RATES)
%plot(CRACK(k).XX,(CRACK(k).PROPAGATION_RATES).*ones(size(CRACK(k).XX)))
%plot(CRACK(k).XX,3.*ones(size(CRACK(k).XX))./3.6)
xl2 = xlabel('Longitude [deg]'); xlim([min(xx),max(xx)]);

tit = title('Cycloid 1');
yl2 = ylabel('Propagation rate [m/s]');
%le2 = legend ('Model','Mean Results','Hoppa 1999');

set(tit,'Fontsize',9,'Interpreter','Latex')
set(xl2,'Fontsize',9,'Interpreter','Latex')
set(yl2,'Fontsize',9,'Interpreter','Latex')
%set(le2,'Fontsize',9,'Interpreter','Latex','Location','NorthWest')

set(gca,'Fontsize',9)
% 
% Cycloid 2
subplot(2,2,3), hold on, grid on, k = 2;
xx = S(k).X; xx = xx(~isnan(xx));
yy = S(k).Y; yy = yy(~isnan(yy));
plot(CRACK(k).XX(1:length(CRACK(k).PROPAGATION_RATES)),CRACK(k).PROPAGATION_RATES)
% plot(CRACK(k).XX,mean(CRACK(k).PROPAGATION_RATES).*ones(size(CRACK(k).XX)))
% plot(CRACK(k).XX,3.*ones(size(CRACK(k).XX)))
xl2 = xlabel('Longitude [deg]'); xlim([min(xx),max(xx)]);

tit = title('Cycloid 2');
yl2 = ylabel('Propagation rate [m/s]');
%le2 = legend ('Model','Mean Results','Hoppa 1999');

set(tit,'Fontsize',9,'Interpreter','Latex')
set(xl2,'Fontsize',9,'Interpreter','Latex')
set(yl2,'Fontsize',9,'Interpreter','Latex')
%set(le2,'Fontsize',9,'Interpreter','Latex','Location','NorthEast')

set(gca,'Fontsize',9)

% % Cycloid 3
% subplot(3,2,5), hold on, grid on, k = 3;
% xx = S(k).X; xx = xx(~isnan(xx));
% yy = S(k).Y; yy = yy(~isnan(yy));
% plot(CRACK(k).XX(1:length(CRACK(k).PROPAGATION_RATES)),CRACK(k).PROPAGATION_RATES)
% % plot(CRACK(k).XX,mean(CRACK(k).PROPAGATION_RATES).*ones(size(CRACK(k).XX)))
% % plot(CRACK(k).XX,3.*ones(size(CRACK(k).XX)))
% xl2 = xlabel('Longitude [deg]'); xlim([min(xx),max(xx)]);
% 
% tit = title('Cycloid 3');
% yl2 = ylabel('Propagation rate [m/s]');
% %le2 = legend ('Model','Mean Results','Hoppa 1999');
% 
% set(tit,'Fontsize',9,'Interpreter','Latex')
% set(xl2,'Fontsize',9,'Interpreter','Latex')
% set(yl2,'Fontsize',9,'Interpreter','Latex')
% %set(le2,'Fontsize',9,'Interpreter','Latex','Location','NorthEast')
% 
% set(gca,'Fontsize',9)

% Delphi Flexus
subplot(2,2,2), hold on, grid on, k = 4;
xx = S(k).X; xx = xx(~isnan(xx));
yy = S(k).Y; yy = yy(~isnan(yy));
plot(CRACK(k).XX(1:length(CRACK(k).PROPAGATION_RATES)),CRACK(k).PROPAGATION_RATES)
% plot(CRACK(k).XX,mean(CRACK(k).PROPAGATION_RATES).*ones(size(CRACK(k).XX)))
% plot(CRACK(k).XX,3.*ones(size(CRACK(k).XX)))
xl2 = xlabel('Longitude [deg]'); xlim([min(xx),max(xx)]);

tit = title('Delphi Flexus');
%yl2 = ylabel('Propagation rate [km/hr]');
%le2 = legend ('Model','Mean Results','Hoppa 1999');

set(tit,'Fontsize',9,'Interpreter','Latex')
set(xl2,'Fontsize',9,'Interpreter','Latex')
%set(yl2,'Fontsize',9,'Interpreter','Latex')
%set(le2,'Fontsize',9,'Interpreter','Latex','Location','NorthWest')

set(gca,'Fontsize',9)

% % Sidon Flexus
% subplot(2,2,4), hold on, grid on, k = 5;
% xx = S(k).X; xx = xx(~isnan(xx));
% yy = S(k).Y; yy = yy(~isnan(yy));
% plot(CRACK(k).XX(1:length(CRACK(k).PROPAGATION_RATES)),CRACK(k).PROPAGATION_RATES)
% % plot(CRACK(k).XX,mean(CRACK(k).PROPAGATION_RATES).*ones(size(CRACK(k).XX)))
% % plot(CRACK(k).XX,3.*ones(size(CRACK(k).XX)))
% xl2 = xlabel('Longitude [deg]'); xlim([min(xx),max(xx)]);
% 
% tit = title('Sidon Flexus');
% %yl2 = ylabel('Propagation rate [km/hr]');
% %le2 = legend ('Model','Mean Results','Hoppa 1999');
% 
% set(tit,'Fontsize',9,'Interpreter','Latex')
% set(xl2,'Fontsize',9,'Interpreter','Latex')
% %set(yl2,'Fontsize',9,'Interpreter','Latex')
% %set(le2,'Fontsize',9,'Interpreter','Latex','Location','NorthEast')
% 
% set(gca,'Fontsize',9)

% Cilicia Flexus
subplot(2,2,4), hold on, grid on, k = 6;
xx = S(k).X; xx = xx(~isnan(xx));
yy = S(k).Y; yy = yy(~isnan(yy));
plot(CRACK(k).XX(1:length(CRACK(k).PROPAGATION_RATES)),CRACK(k).PROPAGATION_RATES)
% plot(CRACK(k).XX,mean(CRACK(k).PROPAGATION_RATES).*ones(size(CRACK(k).XX)))
% plot(CRACK(k).XX,3.*ones(size(CRACK(k).XX)))
xl2 = xlabel('Longitude [deg]'); xlim([min(xx),max(xx)]);

tit = title('Cilicia Flexus');
%yl2 = ylabel('Propagation rate [km/hr]');
%le2 = legend ('Model','Mean Results','Hoppa 1999');

set(tit,'Fontsize',9,'Interpreter','Latex')
set(xl2,'Fontsize',9,'Interpreter','Latex')
%set(yl2,'Fontsize',9,'Interpreter','Latex')
%set(le2,'Fontsize',9,'Interpreter','Latex','Location','NorthEast')

set(gca,'Fontsize',9)

figure, scatter (CRACK(1).XX(1:length(CRACK(1).PROPAGATION_RATES)),CRACK(1).YY(1:length(CRACK(1).PROPAGATION_RATES)),8,...
    CRACK(1).STANDBY,'filled'),colorbar
figure, scatter (CRACK(2).XX(1:length(CRACK(2).PROPAGATION_RATES)),CRACK(2).YY(1:length(CRACK(2).PROPAGATION_RATES)),8,...
    CRACK(2).STANDBY,'filled'),colorbar
figure, scatter (CRACK(3).XX(1:length(CRACK(3).PROPAGATION_RATES)),CRACK(3).YY(1:length(CRACK(3).PROPAGATION_RATES)),8,...
    CRACK(3).STANDBY,'filled'),colorbar
figure, scatter (CRACK(4).XX(1:length(CRACK(4).PROPAGATION_RATES)),CRACK(4).YY(1:length(CRACK(4).PROPAGATION_RATES)),8,...
    CRACK(4).STANDBY,'filled'),colorbar
figure, scatter (CRACK(5).XX(1:length(CRACK(5).PROPAGATION_RATES)),CRACK(5).YY(1:length(CRACK(5).PROPAGATION_RATES)),8,...
    CRACK(5).STANDBY,'filled'),colorbar
figure, scatter (CRACK(6).XX(1:length(CRACK(6).PROPAGATION_RATES)),CRACK(6).YY(1:length(CRACK(6).PROPAGATION_RATES)),8,...
    CRACK(6).STANDBY,'filled'),colorbar

%     figure,
%     loglog(CRACK(k).TIME./(3600*24*30),CRACK(k).PROPAGATION_RATES)
%     grid on, xlabel('log(time [months])'),ylabel('log(propagation rate [m/s])')

% %% MAP PLOT
% 
% figure('Units','inches','Position',[0 0 7.25 3.5],'PaperPositionMode','auto'), hold on
% axesm('MapProjection','eqdcylin','FontName','times','FontSize',10,...
%     'MapLatLimit',[-70 70],'MapLonLimit', [0 360],...
%     'LabelFormat','none',...
%     'MLabelLocation',50,'PLabelLocation',50)
% 
% framem on, mlabel ('south'), plabel on
% [lon,lat] = meshgrid(phi,theta);
% geoshow(lat,lon+180,raster);
% tightmap
% 
% for q = 1:6
%     figure(3),hold on
%     h = scatterm (CRACK(q).YY,CRACK(q).XX + 180,40,CRACK(q).PROPAGATION_RATES,'filled');
%     textm (CRACK(q).YY(1),CRACK(q).XX(1) + 180,num2str(q))
% end
% c = colorbar;
% c.FontName = 'Helvetica';
% c.Location = 'EastOutside';

% figure(10),hold on
% set(figure(10),'units','pixels','position',[0,500,1000,500])
% axesm ('eqdcylin', 'Frame', 'on', 'Grid', 'on','MapLatLim',[-75,75],'MapLonLim',[-360,0]);
% framem on, gridm on, mlabel on, plabel on
% [lon,lat] = meshgrid(phi,theta); 
% geoshow(lat,lon+180,raster);
% geoshow([S(:).Y],[S(:).X]+180,...
%     'DisplayType','Point','Marker','o','Color','red')
% 
% for i= 1:6
%     textm ([S(i).Y(1)],[S(i).X(1)]+180,num2str(i))
% end

% %% Bar plot
% VECT = [max(CRACK(1).PROPAGATION_RATES),CRACK(1).LENGTH/CRACK(1).TIME(end-1);
% max(CRACK(2).PROPAGATION_RATES),CRACK(2).LENGTH/CRACK(2).TIME(end-1);
% max(CRACK(5).PROPAGATION_RATES),CRACK(5).LENGTH/CRACK(5).TIME(end-1);
% max(CRACK(6).PROPAGATION_RATES),CRACK(6).LENGTH/CRACK(6).TIME(end-1)];