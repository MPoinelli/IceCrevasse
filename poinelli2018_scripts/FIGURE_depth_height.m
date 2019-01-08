% Mattia Poinelli
% JPL, April 2017     %%%%% UPDATED JULY 2018 FOR THE PAPER %%%%%%
%
% Depth of superficial crevasse in the location of selected features 
% the scripts loads the shape file Features_smoothed.shp and it is possible
% to select SINGLE crevasses
%
% NORMAL STRESS ROTATED TO the existing lineament

clear all,  clc
%% SELECTION OF PLANET

run Europa_physics.m
run Europa_LEFM_parameters.m

time =3*T_europa/4;

%figure('renderer','zbuffer','Units','inches','Position',[0 0 9 4.5],'PaperPositionMode','auto')
subplot (2,2,4)
%Tit = title ('PeriJove');
Tit = title ('3/4 of Orbital Period');
set(Tit,'Fontsize',8)
set(Tit,'Interpreter','latex')

%% GRID EVALUATION

phi      = [0:2:360];      phi_rad = deg2rad(phi);     % longitudes
theta    = 90 - [0:2:180]; theta_rad = deg2rad(theta); % latitudes

% Features loading
S = shaperead('Features_smoothed.shp');
NAME = struct2cell(S);

%% FEATURE LOOP
for q = [1,2,5,6,12,15,16,21,24] % [1 2 5 6 7 11 12 15 19 21 22 24 26]
    
    % avoid NaN values, xx = longitudes and yy = latitudes
    
    xx = S(q).X; xx = xx(~isnan(xx));
    yy = S(q).Y; yy = yy(~isnan(yy));
    
    % Reuction of size
    xx = xx(1:20:end);
    yy = yy(1:20:end);
    
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
                sin(A),cos(A)];% counter clock-wise rotation (positive angle)
            
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
                sin(-A),cos(-A)];% clock-wise rotation (negative angle)
            
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
        Tangent_stress (i) = T_R(1,1,i);
        Normal_stress (i) = T_R(2,2,i);
        
        %% LEFM loop
        
        % compression -> no fracture a priori
        if Normal_stress(i) <= 0
            depth(i) = 0;
            
            % tension -> LEFM loop
        else
            
            for j = 1 : length(d)
                
                gamma = @(b) b./d(j);
                
                % -----------------------------------------------------------------------------------------------------------
                % Surface Crevasses
                I_overburden = @ (b) (- b + (rho_i -rho_s).* (1 - exp( - C .* b))./ (rho_i * C)) .* G(gamma(b),lambda(j));
                I_waterpress = @ (b) (b - a) .* G(gamma(b),lambda(j));
                
                K_2(j)     = 2*rho_i*g*integral(I_overburden,1,d(j))/sqrt(pi*d(j));
                
%                 % Addition of water pressure filling the crevasse
%                 if d(j) > a
%                     K_3(j) = 2*rho_w*g*integral(I_waterpress,a,d(j))/sqrt(pi*d(j));
%                 else
%                     K_3 (j) = 0;
%                 end
                
                K_1(j) = F(d(j)/H).*Normal_stress(i).*sqrt(pi.*d(j));
                
                K_net(j) = (K_1(j) + K_2(j) )/1000; % [KPa]
%                %-----------------------------------------------------------------------------------------------------------
                
                
%                 % -----------------------------------------------------------------------------------------------------------
%                 % Bottom Crevasses
%                 
%                 sigma_n  = @(z) -rho_i * g .*(H - z) + (rho_i-rho_s)*g.*(1-exp(-C.*(H-z)))./C +rho_w*g.*(H_p-z) + Normal_stress(i);
%                 I_bottom = @(z) 2*G(gamma(z),lambda(j)).*sigma_n(z)./sqrt(pi.*d(j));
%                 
%                 K_net(j) = integral(I_bottom,0,d(j))/1000; % [KPa]
%                 % -----------------------------------------------------------------------------------------------------------
                
            end
            
            location(i) = knnsearch(K_net'- TOU,0); % close value toughness
                        
            % Check for stable solution
            if location(i) == 1
                
                K_net(location(i)) = K_net(location(i) + 1) ;
                location(i) = knnsearch(K_net'- TOU,0); % close value toughness conversion in KPa
                depth(i) = d(location(i));
                
            elseif location(i) == length(K_net)
                error ('Critical depth is higher than the value posed as extreme limit for d')
                return
                
            elseif K_net(location(i) + 1 ) > TOU && K_net(location(i) - 1) < TOU % negative derivative LEFT + RIGHT - for stability
                
                K_net(location(i)) = K_net(location(i) + 1) ;
                location(i) = knnsearch(K_net'- TOU,0); % close value toughness conversion in KPa
                depth(i) = d(location(i));
                
            elseif K_net(location(i) + 1 ) < TOU && K_net(location(i) - 1) < TOU && K_net(location(i)) < TOU
                location (i) = 1;
                depth(i) = d(location(i));
                
            end
            
            depth(i) = d(location(i));
            
        end
    end
    
    if q == 1
    %% Load Raster & plot map
    %[raster] = geotiffread('Europa_raster');
    %[raster] = imread('Europa_raster_fine.tif');
    [raster] = imread('Europa_raster_coarse.jpg');
    
%     axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
%         'MapLatLimit',[-50 -25],'MapLonLimit', [95 200],...
%         'LabelFormat','none',...
%         'MLabelLocation',20,'PLabelLocation',5)

     axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
    'MapLatLimit',[-60 60],'MapLonLimit', [0 360],...
    'LabelFormat','none',...
    'MLabelLocation',50,'PLabelLocation',50)
    
    framem on, mlabel ('south'), plabel on
    [lon,lat] = meshgrid(phi,theta);
    geoshow(lat,lon+180,raster);
    tightmap
    fclose('all')
    %quiverm (yy(1:20:end),xx(1:20:end)+180,Norm_ve(1,1:20:end),Norm_ve(2,1:20:end),'r')
    
    end
    h = scatterm (yy,xx+180,2,depth,'filled');
    %textm(yy(1),xx(1)+180,num2str(q),'white')
    drawnow

    Max_depth (q) = max(depth);
    
    clear xx yy A B T ROT T_R Norm_ve Tang_ve K_net location depth K_1 K_2 K_3 Normal_stress
    
end

max(Max_depth)
% % 
c = colorbar;
c.FontName = 'Times';
c.Location = 'manual'
%c.Location = 'northoutside'
c.Position = [1.5    0.1577    0.0310    0.7100]
c.TickLabels = {'0 m','10 m','20 m','30 m','40 m','50 m','60 m','70 m','80 m','90 m','100 m'};
% c.TickLabels = {'0 m','200 m','400 m','600 m','800 m','1000 m','1200 m'};

% 
 caxis([0 100]);

% set (gca,'Units','Normalized','FontUnits','Points','FontWeight','normal',...
%     'FontSize',11,'FontName','Times')
% set(gcf,'renderer','Painters')
% %%%%%% SLIDES %%%%%%%%%
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% print -dpdf censimento

% % 
% %%%%%% PAPER %%%%%%%%%
% print -depsc2 local_depths_surface

%69 73 86 97
%842 892 1042 1192

% print -depsc2 FIGURE_average_depth
%print(gcf,'-dtiff','ciao.tif')