% Mattia Poinelli
% JPL, April 2017
% 
% Depth of superficial crevasse in the location of selected features 
% the scripts loads the shape file Features_smoothed.shp
% in the time domain
% 
% NORMAL STRESS ROTATED TO the existing lineament
% i.e. at the local scale
%
% KEY OUTPUT:
%
% Max_depth(or height): maximum depth (or height) of the crevasse
% TO MODEL BOTTOM: uncomment 182-189 and comment 162-179
%
% Plot of linement depth at 
% TO PLOT: choose a single value of time and uncomment lines 29-52 and 223-226

%% Load parameters & shape file
run Europa_physics.m
run Europa_LEFM_parameters.m

S = shaperead('Features_smoothed.shp');
NAME = struct2cell(S);
%% TIME STEP of the integration

% input single value for time to plot 
time = linspace (0,T_europa,8);
time = 0;

% %% FIGURE initialization
% 
% phi      = [0:2:360];      phi_rad = deg2rad(phi);     % longitudes
% theta    = 90 - [0:2:180]; theta_rad = deg2rad(theta); % latitudes
% 
% figure('renderer','zbuffer','Units','inches','Position',[0 0 9 4.5],'PaperPositionMode','auto')
% %subplot (2,2,4)
% Tit = title ('PeriJove');
% %Tit = title ('3/4 of Orbital Period');
% set(Tit,'Fontsize',8)
% set(Tit,'Interpreter','latex')
% 
% [raster] = imread('Europa_raster_coarse.jpg');
% 
% axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
%     'MapLatLimit',[-60 60],'MapLonLimit', [0 360],...
%     'LabelFormat','none',...
%     'MLabelLocation',50,'PLabelLocation',50)
% 
% framem on, mlabel ('south'), plabel on
% [lon,lat] = meshgrid(phi,theta);
% geoshow(lat,lon+180,raster);
% tightmap
% fclose('all')
            
%% FEATURE LOOP
for q = 1%[1,2,5,6,12,15,16,21,24] % [1,2,5,6,9,12,15,16]
    feature = q
    % avoid NaN values, xx = longitudes and yy = latitudes
    
    xx = S(q).X; xx = xx(~isnan(xx));
    yy = S(q).Y; yy = yy(~isnan(yy));
    
    % Reuction of size
    xx = xx(1:1:end);
    yy = yy(1:1:end);

        % FOR EVERY TIMESTEP
        for t = 1: length(time)
           time_step = t
           
        % FOR EVERY POINT OF THE FEATURE
        for i = 1 : length(xx)

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

            % Diurnal and secular stress evaluation as Jara-Orue, 2011
            [sigma_theta_day,sigma_phi_day,tau_day] = ...
                diurnal_stress (deg2rad(90 - yy(i)),deg2rad(xx(i)),...
                n,e,w,time(t),epsilon,g,R,mu,eta,h_d,l_d);

            [sigma_theta_sec,sigma_phi_sec,tau_sec] = secular_stress (deg2rad(90 - yy(i)),deg2rad(xx(i)),...
                n,time(t),T_ns,g,R,mu,eta,h_s,l_s);

            % Cauchy Stress Tensor for every point of the Feature
            T(:,:,i) = [sigma_phi_day+sigma_phi_sec,tau_day+tau_sec;
                tau_day+tau_sec,sigma_theta_day+sigma_theta_sec];  % [Pa]

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


%                     % -----------------------------------------------------------------------------------------------------------
%                     % Bottom Crevasses
%                     
%                     sigma_n  = @(z) -rho_i * g .*(H - z) + (rho_i-rho_s)*g.*(1-exp(-C.*(H-z)))./C +rho_w*g.*(H_p-z) + Normal_stress(i);
%                     I_bottom = @(z) 2*G(gamma(z),lambda(j)).*sigma_n(z)./sqrt(pi.*d(j));
%                     
%                     K_net(j) = integral(I_bottom,0,d(j))/1000; % [KPa]
%                     % -----------------------------------------------------------------------------------------------------------

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
        Max_feature_depth_now(t) = max (depth); 
        end

%         %% FIGURE 
%         h = scatterm (yy,xx+180,2,depth,'filled');
%         %textm(yy(1),xx(1)+180,num2str(q),'white')
%         drawnow
        
        Max_depth (q) = max(Max_feature_depth_now);

    clear xx yy A B T ROT T_R Norm_ve Tang_ve K_net location depth K_1 K_2 K_3 Normal_stress Max_feature_depth_now
    
end
