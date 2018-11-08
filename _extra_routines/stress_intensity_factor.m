% Mattia Poinelli
% JPL, Pasadena, March 2017
%
% Application of Linear Elastic Fracture Mechanics to pinpoints Europa in
% order to study stress intensity factor
% Surface & Bottom crevasses possible

clear all, clc, 

%% Load parameters 
run Europa_physics.m
run Europa_LEFM_parameters.m

time = 2*T_europa/4; % [hours]

% pinpoint TARGET: Pywyll Crater
theta =  -15;  % [deg]
phi   =  -270.1; % [deg]

% Stress determination
[sigma_theta_day,sigma_phi_day,tau_day] = ...
    diurnal_stress (deg2rad(90 - theta),deg2rad(phi),...
    n,e,w,time,epsilon,g,R,mu,eta,h_d,l_d);

[sigma_theta_sec,sigma_phi_sec,tau_sec] = secular_stress (deg2rad(90 - theta),deg2rad(phi),...
    n,time,T_ns,g,R,mu,eta,h_s,l_s);

% Tensor for every Feature
T(:,:) = [sigma_theta_day + sigma_theta_sec,tau_day + tau_sec;
    tau_day + tau_sec,sigma_phi_day + sigma_phi_sec];  % [Pa]

% Principal Stresses
sigma_principal(:) = eig(T);

if sigma_principal(2) > 0 % only for tensile status
    
    %% SURFACE& BOTTOM CREVASSES
    % Depth&Height Estimation
    for j = 1 : length(d)
        
        gamma = @(b) b./d(j);
        
        %% Surface Crevasses
        I_overburden = @ (b) (- b + (rho_i -rho_s).* (1 - exp( - C .* b))./ (rho_i * C)) .* G(gamma(b),lambda(j));
        I_waterpress = @ (b) (b - a) .* G(gamma(b),lambda(j));
        
        K_2(j)     = 2*rho_i*g*integral(I_overburden,1,d(j))/sqrt(pi*d(j));
        
% %         % Addition of water pressure filling the crevasse
% %         if d(j) > a
% %             K_3(j) = 2*rho_w*g*integral(I_waterpress,a,d(j))/sqrt(pi*d(j));
% %         else
% %             K_3 (j) = 0;
% %         end
        
        K_1(j) = F(d(j)/H).*sigma_principal(2).*sqrt(pi.*d(j));
        
        K_net_surface(j) = (K_1(j) + K_2(j) )/1000; % [KPa]
        
%         %% Bottom Crevasses
% 
%         sigma_n  = @(z) -rho_i * g .*(H - z) + (rho_i-rho_s)*g.*(1-exp(-C.*(H-z)))./C +rho_w*g.*(H_p-z) + sigma_principal(2);
%         I_bottom = @(z) 2*G(gamma(z),lambda(j)).*sigma_n(z)./sqrt(pi.*d(j));
%         
%         K_net_bottom(j) = integral(I_bottom,0,d(j))/1000;    % [KPa]
        
    end
end

% % LEFM results
location = knnsearch(K_net_surface'- TOU,0); % close value toughness conversion in KPa
            
% Check for stable solution
if location == 1
    
    K_net_surface(location) = K_net_surface(location + 1) ;
    location = knnsearch(K_net_surface'- TOU,0); % close value toughness conversion in KPa
    depth = d(location);
    
elseif location == length(K_net_surface)
    error ('Critical depth is higher than the value posed as extreme limit for d')
    return
    
elseif K_net_surface(location + 1 ) > TOU && K_net_surface(location - 1) < TOU % negative derivative LEFT + RIGHT - for stability
    
    K_net_surface(location) = K_net_surface(location + 1) ;
    location = knnsearch(K_net_surface'- TOU,0); % close value toughness conversion in KPa
    depth = d(location);
    
elseif K_net_surface(location + 1 ) < TOU && K_net_surface(location - 1) < TOU && K_net_surface(location) < TOU
    location  = 1;
    depth = d(location);
    
end

depth = d(location);

%% SURFACE BOTTOM CREVASSES STRESS INTENSITY FACTOR
% figure(1)
% %set(figure(1),'units','pixels','position',[0,500,2200,1250])
% subplot(2,2,4),hold on,grid on
% Tit = title ('3/4 of Orbital Period');
% set(Tit,'Fontsize',15)
% 
% yyaxis left
% plot (K_net_surface,d,'Linewidth',1.5)
% ylim([0 90])
% ylabel('Depth [m]','FontSize',13)
% xlabel('Stress Intensity Factor [ KPa m\^1/2 ]','FontSize',13)
% 
% yyaxis right
% plot (K_net_bottom,d,'Linewidth',1.5)
% ylim([0 1000])
% ylabel('Height [m]','FontSize',13)
% xlim([-200 100])
% % ylabel('h [m]'),xlabel('K [ KPa m\^1/2 ]')

%% ICE THICKNESS INFLUENCE
figure(1),hold on, grid on
set(figure(1),'Units','inches','Position',[0 0 4.5 4.5],'PaperPositionMode','auto')

plot (K_net_surface,d,'Linewidth',1.1)
plot (TOU.*ones(size(0:H)),0:H,'--k','Linewidth',.5)
xlim([-1000, 1000])
ylim([0 60])

% legend('100 m','200 m','500 m', '1000 m', '5000 m', '10000 m','Toughness',...
%     'Location','SouthWest')
legend('1/4 of Orbital Period','1/2 of Orbital Period', 'Toughness','Location','NorthWest')
grid on
ylabel('Depth [m]','FontUnits','Points','Interpreter','Latex',...
    'FontSize',12,'FontName','Times')
xlabel('Stress Intensity Factor [ KPa m$^{1/2}$ ]','FontUnits','Points','Interpreter','Latex',...
    'FontSize',12,'FontName','Times')
set (gca,'Units','Normalized','FontUnits','Points','FontWeight','normal',...
    'FontSize',11,'FontName','Times')
% 
% %%%%%% SLIDES %%%%%%%%%
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% print -dpdf stress_intensity_factor_1
% 
% %%%%%% PAPER %%%%%%%%%
% print -depsc2 stress_intensity_factor_1