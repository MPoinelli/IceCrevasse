% Mattia Poinelli
% JPL, december 2018
%
% Tracking propagation rates

%clear all,  clc, close all
% profile on

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
    (cos(phi).*sin(w + n.*t + atan(Lambda)))) + ... % beginning of secular parts
    C .* alpha_theta_theta (cotheta,2,h_s,l_s) .* ...
    cos(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

sigma_phi = @(t,cotheta,phi) Z .* (-6.*e.*beta_phi_phi(cotheta,0,h_d,l_d).*cos(n.*t+atan(Lambda)) + ...
    e.*beta_phi_phi(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* sin (n.*t + atan(Lambda)) + ...
    3.*cos(2.*phi).*cos(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_phi_phi(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*sin(w + n.*t + atan(Lambda)))) + ...
    C .* alpha_phi_phi (cotheta,2,h_s,l_s) .* ... % beginning of secular parts
    cos(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

tau = @(t,cotheta,phi) Z.* (2.*e.*beta_theta_phi(cotheta,2,h_d,l_d) .* ...
    (4.* cos(2.*phi).* sin(n.*t + atan(Lambda)) - ...
    3.* sin(2.*phi).*cos(n.*t+atan(Lambda))) + ...
    4 .* cos(epsilon).*sin(epsilon).* beta_theta_phi(cotheta,1,h_d,l_d).* ...
    sin(phi).*sin(w + n.*t + atan(Lambda))) - ... % beginning of secular parts
    C .* alpha_theta_phi (cotheta,2,h_s,l_s) .* ...
    sin(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

% Temporal derivative of Stress function
d_sigma_theta = @(t,cotheta,phi) n.* Z .* (6.*e.*beta_theta_theta(cotheta,0,h_d,l_d).*sin(n.*t+atan(Lambda)) + ...
    e.*beta_theta_theta(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* cos (n.*t + atan(Lambda)) - ...
    3.*cos(2.*phi).*sin(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_theta_theta(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*cos(w + n.*t + atan(Lambda)))) - ... % beginning of secular parts
    (4.*pi./T_ns).*C .* alpha_theta_theta (cotheta,2,h_s,l_s) .* ...
    sin(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

d_sigma_phi = @(t,cotheta,phi) n.*Z .* (6.*e.*beta_phi_phi(cotheta,0,h_d,l_d).*sin(n.*t+atan(Lambda)) + ...
    e.*beta_phi_phi(cotheta,0,h_d,l_d).*(4.*sin(2.*phi).* cos (n.*t + atan(Lambda)) - ...
    3.*cos(2.*phi).*sin(n.*t+atan(Lambda))) + ...
    4.*cos(epsilon).*sin(epsilon).* beta_phi_phi(cotheta,1,h_d,l_d) .* ...
    (cos(phi).*cos(w + n.*t + atan(Lambda)))) + ...
    (4.*pi./T_ns).*C .* alpha_phi_phi (cotheta,2,h_s,l_s) .* ... % beginning of secular parts
    sin(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

d_tau = @(t,cotheta,phi) n.*Z.* (2.*e.*beta_theta_phi(cotheta,2,h_d,l_d) .* ...
    (4.* cos(2.*phi).* cos(n.*t + atan(Lambda)) + ...
    3.* sin(2.*phi).*sin(n.*t+atan(Lambda))) + ...
    4 .* cos(epsilon).*sin(epsilon).* beta_theta_phi(cotheta,1,h_d,l_d).* ...
    sin(phi).*cos(w + n.*t + atan(Lambda))) - ... % beginning of secular parts
    (4.*pi./T_ns).*C .* alpha_theta_phi (cotheta,2,h_s,l_s) .* ...
    cos(2.*phi + 4 .* pi .* t ./T_ns + atan(Delta));

%% MAPPING SETTINGS

% Features loading
S = shaperead('Features_smoothed.shp');

%depth (1,:) =  max(max(depth)); depth (end,: ) =  max(max(depth));
load('Critical_depths_file.mat')
depth = depth(:,:,1:end-2);
Mean_depth = mean(depth,3);

% map
% figure('Units','inches','Position',[0 0 9 4.5],'PaperPositionMode','auto')
% 
% axesm ('eqdcylin');
% axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
%     'MapLatLimit',[-70 70],'MapLonLimit', [10 350],...
%     'LabelFormat','none',...
%     'MLabelLocation',30,'PLabelLocation',30)
% tightmap
% 
% mlabel ('south'), plabel on
% [lon,lat] = meshgrid(phi,theta);
% geoshow(lat,lon+180,raster);

% generation of coarse grid
% useful only when dealing with mean depths

phi      = [0:20:360]; phi_rad = deg2rad(phi);          % longitudes
theta    = 90 - [0:20:180]; theta_rad = deg2rad(theta); % latitudes
%  contourm (theta,phi,Mean_depth,'LevelStep',3,'ShowText','off','Linewidth',1.3);
% % 
%  geoshow([S(1:end).Y],[S(1:end).X]+180,...
%      'Color','black','Linewidth',2)
for q = 1:26 % 1:26% [1,2,5,6,12,15,16,21,24] %  selection of a the feature from the shape file
    
    time = 0;
    TIME = [];
    
    % avoid NaN values, xx = longitudes and yy = latitudes
    xx = S(q).X; xx = xx(~isnan(xx));
    yy = S(q).Y; yy = yy(~isnan(yy));
    END = length(xx);
    Perc_completion (q) = 100;
    
    CRACK(q).XX = xx;
    CRACK(q).YY = yy;
    
    % loop on the feature
    for i = 1:END
        
        START = time;
        check = 1;
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
        end
                
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
            Normal_stress   (i) = T_R(2,2,i);
            d_Normal_stress (i) = d_T_R(2,2,i);
            Shear_stress    (i) = T_R(1,2,i);
            d_Shear_stress  (i) = d_T_R(1,2,i);
            
            % Normal displacement rate
            epsilon_dot     (i) = d_STRAIN_R(2,2,i);
            epsilon(i) = STRAIN_R(2,2,i);
            
            % Normal stress for the rest of the feature at this Time
            % Multiplication of Matrices
            N_S_feat = ROT(2,1,:).^2.*T(1,1,:)+ROT(2,1,:).*ROT(2,2,:).*(T(2,1,:)+T(1,2,:))+...
                ROT(2,2,:).^2.*T(2,2,:);
               
            if i ~= 1              
                for j = i-1:-1:1 % check which part of the crack is active
                    
                    if N_S_feat(j) > 0
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
            
            % Propagation rate from the derivative of displacement control
            V_rate(i) =  (E/4*opening_rate(i) - d_Normal_stress(i)*a*V(a,b))/...
                (Normal_stress(i)*V(a,b) + Normal_stress(i)*a*dV(a,b));
            
            %% PROPAGATION CHECK
                            
            if K_I(i) > TOU*1000 && V_rate(i) > 0 
                % propagation
                check = -1;
                
            elseif time - START >= 20* T_europa
                STRING = strcat('Feature formed at :',num2str(floor(i/END*100)),'%');
                %disp(STRING)
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
                time = time + 1800;
                clear a 
            end
            
        end
        
        % time step for finalising the shaped segment
        time_step(i) = angular_distance(i)*R/2/V_rate(i);
        
        time = time + time_step(i);
        TIME = [TIME,time];
       
        if i == END % End of the feature
            CRACK(q).PROPAGATION_RATES = V_rate;
            CRACK(q).OPENING_RATE = opening_rate;
            CRACK(q).WIDTH = delta;
            CRACK(q).TIME = TIME;
            CRACK(q).EPSILON_DOT = epsilon_dot;
            XX = xx;
            YY = yy;
            
            % Problem with Phoenix linea
            if q ==24
               time = TIME(end)-TIME(1);
            end
            
            LEGEND{q} = num2str(S(q).Name);
            Cycles_to_build(q) = time/T_europa;
            Max_opening_width(q) = max(delta);
            Max_opening_rate (q) = max(opening_rate); % m/s
            Max_prop_rate(q) = max(V_rate);
            Mean_prop_rate (q) = a/time;%*3600/1000; 
                       
            clear a delta opening_rate time_step V_rate Normal_stress epsilon_dot epsilon 
            clear SIGMA_COTHETA SIGMA_PHI TAU STRAIN
            clear D_SIGMA_COTHETA D_SIGMA_PHI D_TAU d_STRAIN
            clear xx yy
            break
        end
        
    end
    
end
% 

k = 1;
xx = S(k).X; xx = xx(~isnan(xx));
yy = S(k).Y; yy = yy(~isnan(yy));
    
figure,hold on,title('zPropagation rate')
plot(xx,CRACK(k).PROPAGATION_RATES)
