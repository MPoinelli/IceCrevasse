% Mattia Poinelli
% JPL, may 2017
%
% Fourier Analyis for the depth signal
% Implementation of FFT for a SINGLE crevasse depth signal in time 

clear all,  clc, close all

% profile on
%% LOADING DATA

% Parameters
run Europa_physics.m
run Europa_LEFM_parameters.m

time = 1:1800:1*T_europa;  % Timestep half an hour

%% GRID EVALUATION
phi      = [0:30:360]; phi_rad = deg2rad(phi);          % longitudes
theta    = 90 - [30:30:150]; theta_rad = deg2rad(theta); % latitudes

% Features loading
S = shaperead('Features_smoothed.shp');

% 1 Agenor
% 8 Mehen Linea
% 9 Agave Linea
% 15 Tormsdale Linea

xx = S(15).X(1:50:end); xx = xx(~isnan(xx));
yy = S(15).Y(1:50:end); yy = yy(~isnan(yy));

for t = 1 :length(time)
    t
    %% Loop for every grid point
    
    % FOR EVERY POINT OF THE FEATURE
    for i = 1 :length(xx)
        
        % Diurnal and secular stress evaluation as Jara-Orue, 2011
        [sigma_theta_day,sigma_phi_day,tau_day] = ...
            diurnal_stress (deg2rad(90 - yy(i)),deg2rad(xx(i)),...
            n,e,w,time(t),epsilon,g,R,mu,eta,h_d,l_d);
        
        [sigma_theta_sec,sigma_phi_sec,tau_sec] = secular_stress (deg2rad(90 - yy(i)),deg2rad(xx(i)),...
            n,time(t),T_ns,g,R,mu,eta,h_s,l_s);
        
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
            
        elseif yy(i+1) > yy(i)% 'ascending' feature 
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
        
        % Normal and Shear Stress
        Shear_stress (i,t) = T_R(1,2,i);
        Normal_stress (i,t) = T_R(2,2,i);
        
        
        %% LEFM analysis
        
        % compression -> no fracture a priori
        if Normal_stress(i,t) <= 0
            
            depth(i,t) = 0;
            
            % tension -> inititalize LEFM loop
        else
            
            for j = 1 : length(d)
                
                gamma = @(b) b./d(j);
                
                %% Surface Crevasses
                I_overburden = @ (b) (- b + (rho_i -rho_s).* (1 - exp( - C .* b))./ (rho_i * C)) .* G(gamma(b),lambda(j));
                I_waterpress = @ (b) (b - a) .* G(gamma(b),lambda(j));
                
                K_2(j)     = 2*rho_i*g*integral(I_overburden,1,d(j))/sqrt(pi*d(j));
                
%                 % Addition of water pressure filling the crevasse
%                 if d(j) > a
%                     K_3(j) = 2*rho_w*g*integral(I_waterpress,a,d(j))/sqrt(pi*d(j));
%                 else
%                     K_3 (j) = 0;
%                 end
                
                K_1(j) = F(d(j)/H).*Normal_stress(i,t).*sqrt(pi.*d(j));
                
                K_net(j) = (K_1(j) + K_2(j) )/1000;
                
                %                 %% Bottom Crevasses
                %
                %                 % Bottom Crevasses
                %                 sigma_n  = @(z) -rho_i * g .*(H - z) + (rho_i-rho_s)*g.*(1-exp(-C.*(H-z)))./C +rho_w*g.*(H_p-z) + Normal_stress(i,t);
                %                 I_bottom = @(z) 2*G(gamma(z),lambda(j)).*sigma_n(z)./sqrt(pi.*d(j));
                %
                %                 K_net(j) = integral(I_bottom,0,d(j))/1000; % [KPa]
                
            end
            
            location(i,t) = knnsearch(K_net'- TOU,0); % close value toughness
            
            % Check for stable solution
            if location(i,t) == 1
                
                K_net(location(i,t)) = K_net(location(i,t) + 1) ;
                location(i,t) = knnsearch(K_net'- TOU,0); % close value toughness conversion in KPa
                depth(i,t) = d(location(i,t));
                
            elseif location(i,t) == length(K_net)
                error ('Critical depth is higher than the value posed as extreme limit for d')
                return
                
            elseif K_net(location(i,t) + 1 ) > TOU && K_net(location(i,t) - 1) < TOU % negative derivative LEFT + RIGHT - for stability
                
                K_net(location(i,t)) = K_net(location(i,t) + 1) ;
                location(i,t) = knnsearch(K_net'- TOU,0); % close value toughness conversion in KPa
                depth(i,t) = d(location(i,t));
                
            elseif K_net(location(i,t) + 1 ) < TOU && K_net(location(i,t) - 1) < TOU && K_net(location(i,t)) < TOU
                location (i,t) = 1;
                depth(i,t) = d(location(i,t));
                
            end
            
            depth(i,t) = d(location(i,t));
            
%             % CHECK SOLUTION
%             figure(2),hold on
%             plot (K_net,d)  
%             plot (K_net(location(i,q)),depth(i,q),'o')
%             drawnow
            
        end
        
    end
end

Max_val     = max(depth);
Mean_depth  = mean(depth,1);
Min_val     = min(depth);

Mean_N_stress = mean(Normal_stress);
Mean_T_stress = mean(Shear_stress);

figure('Units','inches','Position',[0 0 6 4.5],'PaperPositionMode','auto')
hold on, grid on

%subplot(1,2,1),hold on,grid on
% T1 = title ('Mean Depth and Normal Stress for Tormsdale Linea');
% set(T1,'Fontsize',15)
yyaxis left
plot (time./(24*3600), Mean_depth,'Linewidth',1.4)
plot (time./(24*3600), Max_val,'Linewidth',.9)
%plot (time./3600, Min_val)
X1 = xlabel('Days after PeriJove');
Y1 = ylabel('Depth [m]');
xlim([0,3.55])
yyaxis right
plot (time./(24*3600),Mean_N_stress./1000,'Linewidth',1.4)
%plot (time./(24*3600),Mean_T_stress./1000,'Linewidth',1.4)
Y2 = ylabel(' Stress [KPa]','FontSize',13);

L = legend('Mean Depth', 'Maximum Depth','Normal Stress');

set(X1,'Fontsize',20);
set(X1,'Interpreter','latex');
set(Y1,'Fontsize',20);
set(Y1,'Interpreter','latex');
set(Y2,'Fontsize',20);
set(Y2,'Interpreter','latex');
set(L ,'Fontsize',15);
set(L ,'Interpreter','latex');
% set(T1,'Fontsize',20);
% set(T1,'Interpreter','latex');
% 
% % Load Raster & plot map
% [raster] = geotiffread('Europa_raster_fine');
% S = shaperead('Features_smoothed.shp');
% 
% %depth (1,:) =  max(max(depth)); depth (end,: ) =  max(max(depth));
% 
% figure(1),hold on
% %set(figure(1),'units','pixels','position',[0,500,2200,1250])
% %subplot(2,2,4)
% Tit = title ('Averaged Critical depth');
% set(Tit,'Fontsize',15)
% %set(Tit,'Interpreter','latex')
% 
% axesm ('eqdcylin','MapLatLim',[-60,60],'MapLonLim',[-360,0]);
% framem on, mlabel ('south'), plabel on
% [lon,lat] = meshgrid(phi,theta);
% 
% contourfm (theta,phi,Mean_depth,'Linestyle','none');c = colorbar;
% c.Label.String = 'Average Critical Depth [m]';c.FontSize = 12;
% c.Limits = [0 max(max(Mean_depth))];
% c.FontName = 'latex';
% 
% geoshow([S(:).Y],[S(:).X]+180,...
%     'Color','white','Linewidth',2)

%save('Critical_depths_file','depth')

%
% figure (2), hold on
% for o = 1 : length(time),
%     D1(o) = depth(5,10,o);
%     D2(o) = depth(10,10,o)
%
% end
%
% plot(time,D1)
% plot(time,D2)

% %% FFT Analysis
% Ts = 1800;
% Fs = 1/Ts;
% 
% L = length(time);
% Y = fft(Mean_depth);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% 
% subplot(1,2,2),hold on,grid on
% T2 = title ('Fourier Analysis for Tormsdale Linea');
% set(T2,'Fontsize',20)
% plot(f.*1000000,P1,'Linewidth',1.4)
% X1 = xlabel('Frequency ($\mu$Hz)');
% Y1 = ylabel('Power Spectrum [m]');
% 
% set(X1,'Fontsize',20);
% set(X1,'Interpreter','latex');
% set(Y1,'Fontsize',20);
% set(Y1,'Interpreter','latex');
% set(T2,'Fontsize',20);
% set(T2,'Interpreter','latex');

% figure, hold on
% plot (time,depth(3,:))
% plot (time,depth(1,:))
% plot (time,depth(5,:))
% plot (time,depth(7,:))
% profile viewer