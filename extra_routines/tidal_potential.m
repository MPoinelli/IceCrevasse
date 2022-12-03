% Non_zero eccentricyt & NSR
clear all , close all, clc

theta = (1:1:179)*pi/180;
lambda = (0:1:359)*pi/180;

n = deg2rad(101.37472)/(60*60*24); % rad/s
% t = [0,0.5*3.85*60*60];         % sec
t = [38357.6,191788]
r = 1562000;   % m
e = 0.0094;

phi_1 = zeros(length(theta),length(lambda)); 
phi_2 = zeros(length(theta),length(lambda)); 

for i = 1: length(theta)
    
    P_2 = norm_ass_Legendre_function (2,cos(theta(i)));

    
    for j= 1:length(lambda)
    
        phi_1(i,j) = ((-0.5*P_2(2,end)+(1/4)*P_2(2,2)*cos(2*lambda(j))) + ...
            ((-3/2)*P_2(2,end)*cos(n*t(1))+(1/2)*P_2(2,2)*(3*cos(2*lambda(j))*cos(n*t(1))) + ...
            4*sin(2*lambda(j))*sin(n*t(1))));
        
        phi_2(i,j) = (n*r)^2 * ( (-0.5*P_2(2,end)+(1/4)*P_2(2,2)*cos(2*lambda(j))) + ...
            e*((-3/2)*P_2(2,end)*cos(n*t(2))+(1/2)*P_2(2,2)*(3*cos(2*lambda(j))*cos(n*t(2))) + ...
            4*sin(2*lambda(j))*sin(n*t(2))));
    
    end
    
    clear P_2
    
end

cmin = min(min(phi_1)); cmax = max(max(phi_1)); % [Gal* m]
cmin = min(min(phi_2)); cmax = max(max(phi_2)); % [Gal*m]
diff = phi_2 - phi_1;
cmin_diff = min(min(diff)); cmax_diff = max(max(diff)); % [Gal*m]

plot_map_gen(lambda, theta, phi_1 , cmin, cmax);
colorbar;

plot_map_gen(lambda, theta, phi_2 , cmin, cmax);
colorbar;

plot_map_gen(lambda, theta, diff , cmin_diff, cmax_diff);
colorbar;
