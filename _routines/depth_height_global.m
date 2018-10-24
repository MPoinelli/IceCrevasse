% Mattia Poinelli
% JPL, April 2017
% 
% Depth of superficial crevasse at the global scale
% the scripts also loads the shape file Features_smoothed.shp
% 
% Stress source is the most tensile principal stress of Jara Orue &
% Vermeersen 2011
%
% KEY OUTPUT:
% 
% depth(i,q,t) : matrix of critical depth for a single location (i,q) at
% the time t
% TO CALCULATE MEAN DEPTH: add several values in the time vector
% 
% Mean_depth(i,q): average depth
%
% Global Map of critical depth (or height)


clear all,  clc, close all

%% Load parameters & shape file
run ../_physical_parameters/Europa_physics.m
run ../_physical_parameters/Europa_LEFM_parameters.m

S = shaperead('../_target_features/Features_smoothed.shp');
NAME = struct2cell(S);

%% TIME STEP and GRID
%time = 1:5*3600:1*T_europa;  % Timestep 60 mins
time = [1,3600];

phi      = [0:20:360]; phi_rad = deg2rad(phi);          % longitudes
theta    = 90 - [0:20:180]; theta_rad = deg2rad(theta); % latitudes

for t = 1 :length(time)
    t
    %% Loop for every grid point
    for q = 1 : length(phi)
        
        for i = 1 : length(theta)
            
             % Diurnal and secular stress evaluation as Jara-Orue, 2011
            [sigma_theta_day,sigma_phi_day,tau_day] = ...
                diurnal_stress (deg2rad(90 - theta(i)),deg2rad(phi(q)),...
                n,e,w,time(t),epsilon,g,R,mu,eta,h_d,l_d);
            [sigma_theta_sec,sigma_phi_sec,tau_sec] = secular_stress (deg2rad(90 - theta(i)),deg2rad(phi(q)),...
                n,time(t),T_ns,g,R,mu,eta,h_s,l_s);
            
            % Cauchy Stress Tensor for every point of the Feature
            T = [sigma_phi_day+sigma_phi_sec,tau_day+tau_sec;
                tau_day+tau_sec,sigma_theta_day+sigma_theta_sec];  % [Pa]

            % Eigenvalues
            [~,D1] = eig (T);
            E  = [D1(1,1);D1(2,2)];
            
            % Most-Tensile Principal Stress is the E(2)
            Principal_stress_1 (i,q,t) = E(1);
            Principal_stress_2 (i,q,t) = E(2);
            
            clear T E D
            
            %% LEFM analysis
            
            % compression -> no fracture a priori
            if Principal_stress_2(i,q,t) <= 0
                
                depth(i,q,t) = 0;
                
                % tension -> inititalize LEFM loop
            else
                
                for j = 1 : length(d)
                    
                    gamma = @(b) b./d(j);
                    
                                        %% Surface Crevasses
                                        I_overburden = @ (b) (- b + (rho_i -rho_s).* (1 - exp( - C .* b))./ (rho_i * C)) .* G(gamma(b),lambda(j));
%                     
                                        K_2(j)     = 2*rho_i*g*integral(I_overburden,1,d(j))/sqrt(pi*d(j));
%                     
%                                         % Addition of water pressure filling the crevasse
%                                         I_waterpress = @ (b) (b - a) .* G(gamma(b),lambda(j));
%                                         if d(j) > a
%                                             K_3(j) = 2*rho_w*g*integral(I_waterpress,a,d(j))/sqrt(pi*d(j));
%                                         else
%                                             K_3 (j) = 0;
%                                         end
                    
                                        K_1(j) = F(d(j)/H).*Principal_stress_2(i,q,t).*sqrt(pi.*d(j));
                    
                                        K_net(j) = (K_1(j) + K_2(j))/1000;
                    
%                     %% Bottom Crevasses
%                     
%                     % Bottom Crevasses
%                     sigma_n  = @(z) -rho_i * g .*(H - z) + (rho_i-rho_s)*g.*(1-exp(-C.*(H-z)))./C +rho_w*g.*(H_p-z) + Principal_stress_2(i,q,t);
%                     I_bottom = @(z) 2*G(gamma(z),lambda(j)).*sigma_n(z)./sqrt(pi.*d(j));
%                     
%                     K_net(j) = integral(I_bottom,0,d(j))/1000; % [KPa]
                    
                end
                
                location(i,q,t) = knnsearch(K_net'- TOU,0); % close value toughness
                
                % Check for stable solution
                if location(i,q,t) == 1
                    
                    K_net(location(i,q,t)) = K_net(location(i,q,t) + 1) ;
                    location(i,q,t) = knnsearch(K_net'- TOU,0); % close value toughness conversion in KPa
                    depth(i,q,t) = d(location(i,q,t));
                    
                elseif location(i,q,t) == length(K_net)
                    error ('Critical depth is higher than the value posed as extreme limit for d')
                    return
                    
                elseif K_net(location(i,q,t) + 1 ) > TOU && K_net(location(i,q,t) - 1) < TOU % negative derivative LEFT + RIGHT - for stability
                    
                    K_net(location(i,q,t)) = K_net(location(i,q,t) + 1) ;
                    location(i,q,t) = knnsearch(K_net'- TOU,0); % close value toughness conversion in KPa
                    depth(i,q,t) = d(location(i,q,t));
                    
                elseif K_net(location(i,q,t) + 1 ) < TOU && K_net(location(i,q,t) - 1) < TOU && K_net(location(i,q,t)) < TOU
                    location (i,q,t) = 1;
                    depth(i,q,t) = d(location(i,q,t));
                    
                end
                
                depth(i,q,t) = d(location(i,q,t));
                
            end
            
        end
    end
end

Mean_depth = mean(depth,3);

%% FIGURE
%depth (1,:) =  max(max(depth)); depth (end,: ) =  max(max(depth));
figure('Units','inches','Position',[0 0 9 4.5],'PaperPositionMode','auto')
%Tit = title ('Average Critical Depth');
% set(Tit,'Fontsize',8)
% set(Tit,'Interpreter','latex'))
axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
    'MapLatLimit',[-90 90],'MapLonLimit', [0 360],...
    'LabelFormat','none',...
    'MLabelLocation',30,'PLabelLocation',30)

framem on, mlabel ('south'), plabel on
[lon,lat] = meshgrid(phi,theta);

contourfm (theta,phi,Mean_depth,'Linestyle','none');

c = colorbar;
c.FontName = 'Times';
%c.Location = 'manual'
%c.Position = [0.9165    0.1577    0.0310    0.7100]
%c.Label.String = 'Average Critical Depth [m]';c.FontSize = 12;
c.Limits = [0 max(max(Mean_depth))];
%c.TickLabels = {'0 m','10 m','20 m','30 m','40 m','50 m','60 m','70 m'};
tightmap, gridm on
geoshow([S([1,2,5,6,12,15,16,21,24]).Y],[S([1,2,5,6,12,15,16,21,24]).X]+180,...
    'Color','white','Linewidth',2)