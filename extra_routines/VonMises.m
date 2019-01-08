% Mattia Poinelli
% JPL, April 2017
%
% Von Mises criterion applied to Europa using MOST-Tensile principal stress
% as opening source.
% Visualisation of failure areas
% Non rotated tensor

clear all,  clc

%% Load parameters & grid evaluation
run Europa_physics.m

time = 3*T_europa/4;        % [hours]
% ------------------------------------------------------------------------------
phi      = [0:10:360];      phi_rad = deg2rad(phi);     % longitudes
theta    = 90 - [0:10:180]; theta_rad = deg2rad(theta); % latitudes

% Critial Failure Parametres
TOU = 100;    % [KPa m1/2] ice toughness
TS  = 100000; % [Pa]       ice tensile strength 

%% Loop for every grid point
for q = 1 : length(phi)
    
    for i = 1 : length(theta)
        
        %% Stress determination
        [sigma_theta_day,sigma_phi_day,tau_day] = ...
            diurnal_stress (deg2rad(90 - theta(i)),deg2rad(phi(q)),...
            n,e,w,time,epsilon,g,R,mu,eta,h_d,l_d);
        [sigma_theta_sec,sigma_phi_sec,tau_sec] = secular_stress (deg2rad(90 - theta(i)),deg2rad(phi(q)),...
            n,time,T_ns,g,R,mu,eta,h_s,l_s);
        
        % Tensor for every point of the Feature
        T = [sigma_theta_day+sigma_theta_sec,tau_day+tau_sec;
            tau_day+tau_sec,sigma_phi_day+sigma_phi_sec];  % [Pa]
        
        [~,D] = eig (T);
        E  = [D(1,1);D(2,2)];
        
        % Most-Tensile Principal Stress
        Principal_stress_1 (i,q) = E(1);
        Principal_stress_2 (i,q) = E(2);
        
        clear T E D
     
        % Von Mises Results
        J_2(i,q) = (1/3)*(Principal_stress_2(i,q)+Principal_stress_1(i,q)).^2 - ...
            Principal_stress_2(i,q).*Principal_stress_1(i,q);
        
        VM(i,q) = (sqrt(3*J_2(i,q))/TS); % percentage of failure

    end
    
end


%% Load Raster & plot map
[raster] = imread('Europa_raster_coarse.jpg');
S = shaperead('Features_smoothed.shp');

VM (1,:) =  max(max(VM)); VM (end,: ) =  max(max(VM));

figure(1),hold on
%set(figure(1),'units','pixels','position',[0,500,2200,1250])
%subplot(2,2,2)
%title ('Linear Elastic Fracture Mechanics')
Tit = title ('PeriJove');
set(Tit,'Fontsize',15)
%set(Tit,'Interpreter','latex')
colormap (hot)
[lon,lat] = meshgrid(phi,theta);

% title ('Von Mises Criterion')
axesm ('eqdcylin','MapLatLim',[-90,90],'MapLonLim',[-360,0]);
framem on, mlabel ('south'), plabel on
% geoshow(lat,lon,raster);
contourfm (theta,phi,VM,'Linestyle','none');c = colorbar; 
c.Label.String = 'Percentage of Failure';c.FontSize = 12;

c.FontName = 'latex';
%.Limits = [0 100]; 
% c.TickLabels = {'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};

geoshow([S(:).Y],[S(:).X]+180,...
    'Color','cyan','Linewidth',2)