% Mattia Poinelli
% JPL, April 2017

% Visualisation of localised rotation on map and on graphs for
% SINGLE crevasses of the moon. NO LEFM implemented
% Rotated tensor

clear all,  clc, close all

%profile on

%% SELECTION OF PLANET
run Europa_physics.m
run Europa_LEFM_parameters.m

time = 3*T_europa/4;        % [hours]
figure('Units','inches','Position',[0 0 9 4.5],'PaperPositionMode','auto')
axesm ('eqdcylin');
axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
    'MapLatLimit',[-70 70],'MapLonLimit', [0 360],...
    'LabelFormat','none',...
    'MLabelLocation',30,'PLabelLocation',30)

%% GRID EVALUATION

phi      = [0:1:360];      phi_rad = deg2rad(phi);     % longitudes
theta    = 90 - [0:1:180]; theta_rad = deg2rad(theta); % latitudes

% Features loading
S = shaperead('Features_smoothed.shp');

for q = [1,2,5,6,9,12,16] % 1:length(time)
    q
    % avoid NaN values, xx = longitudes and yy = latitudes
    
    xx = S(q).X; xx = xx(~isnan(xx));
    yy = S(q).Y; yy = yy(~isnan(yy));
    
    % FOR EVERY POINT OF THE FEATURE
    for i = 1 : length(xx)
        
        % Diurnal and secular stress evaluation as Jara-Orue, 2011
        [sigma_theta_day,sigma_phi_day,tau_day] = ...
            diurnal_stress (deg2rad(90 - yy(i)),deg2rad(xx(i)),...
            n,e,w,time,epsilon,g,R,mu,eta,h_d,l_d);
        
        [sigma_theta_sec,sigma_phi_sec,tau_sec] = secular_stress (deg2rad(90 - yy(i)),deg2rad(xx(i)),...
            n,time,T_ns,g,R,mu,eta,h_s,l_s);
        
        % Cauchy Stress Tensor for every point of the Feature
        T(:,:,i) = [sigma_phi_day+sigma_phi_sec,tau_day+tau_sec;
            tau_day+tau_sec,sigma_theta_day+sigma_theta_sec];  % [Pa]                      
        
        %% Rotation to local coordinates
        if i == length(xx)
            % last node of the grid has the same charateristic of the
            % previous node to avoid Matlab unconsistencies
            
            ROT(:,:,i) =  ROT(:,:,i-1);
            
            Tang_ve(1,i) = Tang_ve(1,i-1);
            Tang_ve(2,i) = Tang_ve(2,i-1);
            Norm_ve(1,i) = Norm_ve(1,i-1);
            Norm_ve(2,i) = Norm_ve(2,i-1);
            
        elseif yy(i+1) > yy(i) % 'ascending' feature 
            % evaluation of useful angles a,b
            b_ang(i) = deg2rad(xx(i+1) - xx(i));
            a_ang(i) = deg2rad(yy(i+1) - yy(i));
            
            % Side_Angle_Side problem, A is the rotation angle, and the
            % last output c gives the length of the segment
            [A,~,~] = Side_Angle_Side (a_ang(i),b_ang(i));
            
            ROT(:,:,i) = [cos(A),-sin(A);
                sin(A),cos(A)]; % counter clock-wise rotation (positive angle)
            
            % Tangent vector
            Tang_ve(1,i) = ROT(1,1,i);
            Tang_ve(2,i) = ROT(2,1,i);
            
            % Normal vector
            Norm_ve(1,i) = ROT(1,2,i);
            Norm_ve(2,i) = ROT(2,2,i);
   
        elseif yy(i+1) < yy(i)% 'descending' feature 
            % evaluation of useful angles a,b
            b_ang(i) = deg2rad(xx(i+1) - xx(i));
            a_ang(i) = abs(deg2rad(yy(i+1) - yy(i)));
            
            % Side_Angle_Side problem, A is the rotation angle, and the
            % last output c gives the length of the segment
            [A,~,~] = Side_Angle_Side (a_ang(i),b_ang(i));
            
            ROT(:,:,i) = [cos(-A),-sin(-A);
                sin(-A),cos(-A)]; % clock-wise rotation (negative angle)
            
            % Tangent vector
            Tang_ve(1,i) = ROT(1,1,i);
            Tang_ve(2,i) = ROT(2,1,i);
            
            % Normal vector
            Norm_ve(1,i) = ROT(1,2,i);
            Norm_ve(2,i) = ROT(2,2,i);
   
        end
        
        % Stress tensor rotated to LOCAL coordinates
        T_R(:,:,i) = ROT(:,:,i)*T(:,:,i)*ROT(:,:,i)';
        
        % Tangent, Normal and Shear Stress
        Tangen_stress (i,q) = T_R(1,1,i);
        Normal_stress (i,q) = T_R(2,2,i);
        Shear_stress  (i,q) = T_R(1,2,i);
    end
    
    %% Load Raster & plot map
  %  [raster] = geotiffread('Europa_raster_coarse');
   

    [raster] = imread('Europa_raster_coarse.jpg');
   

%     Tit = title ('Directions of Stresses');
%     set(Tit,'Fontsize',15)
%     axesm ('eqdcylin','MapLatLim',[-60,60],'MapLonLim',[-360,0]);
%     framem on, mlabel ('south'), plabel on
    [lon,lat] = meshgrid(phi,theta);
    if q ==1 
    geoshow(lat,lon+180,raster);
    end
    plotm (yy,xx+180,'w','LineWidth',1.5)
    
    quiverm (yy(250),xx(250)+180,7 * Norm_ve(2,250),7 * Norm_ve(1,250),'r')
    quiverm (yy(250),xx(250)+180,7 * Tang_ve(2,250),7 * Tang_ve(1,250),'b')
    
    quiverm (yy(100),xx(100)+180,7 * Norm_ve(2,100),7 * Norm_ve(1,100),'r')
    quiverm (yy(100),xx(100)+180,7 * Tang_ve(2,100),7 * Tang_ve(1,100),'b')
    
    drawnow
    %clear xx yy A B T ROT T_R Norm_ve Tang_ve

end
% % 
% figure (2), hold on
% Tit = title ('Shear and normal stress')
% plot (time./(3600),Normal_stress./1000)
% plot (time/(3600),Shear_stress./1000)
% legend ('Normal Stress','Shear Stress')
% % profile viewer
tightmap
framem on, mlabel ('south'), plabel on
set (gca,'Units','Normalized','FontUnits','Points','FontWeight','normal',...
    'FontSize',11,'FontName','Times')
% %%%%%% SLIDES %%%%%%%%%
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% print('-dpdf','-r100')

% 
% %%%%%% PAPER %%%%%%%%%
% print -depsc2 rotation_stress
%save2pdf ('Rotatedtensor')