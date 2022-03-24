% Europa parameters useful for a general LEFM analysis
% from Van Der Veen, 1998

%% Ice Thickness & fracture Depth

% surface crevasses
global H d

H = 5000;  % [m]
d = 5:5:150; 

% %bottom crevasses
% global H d a 
% H = 5000;  % [m]
% d = 5:20:1500; 
% a = H;     % water penetration

%% Physical Properties of the Ice
global rho_i rho_s rho_w C g

rho_i = 917;  % [kg/m3]
rho_s = 850;  % [kg/m3]
rho_w = 1000; % [kg/m3]
C     = 0.02; % [~]
g     = 1.35; % [m/s2]

% Critial Failure Parametres
global TOU TS

TOU = 100;    % [KPa m1/2] ice toughness
TS  = 100000; % [Pa]       ice tensile strength

% Average Density and Piezometric Head
global rho_a H_p

rho_a = rho_i - (rho_i - rho_s)*(1 - exp(-C*H))/(C*H); % [kg/m3] average density
H_p = rho_a*H/rho_w; % [m]

%% LEFM tools, from Tada 2000
global lambda F G

lambda = d./H;

F = @(x) 1.12 - 0.23 .* x + 10.55 .* x .^2 -21.72 .* x .^3 + ...
    30.39 .* x .^ 4;
G = @(gamma,lambda) (3.52.*(1 - gamma)) ./ (1 - lambda).^(3/2) - ...
    (4.35 - 5.28.*gamma) ./ (1 - lambda).^(1/2) + ( (1.3 - 0.3.*gamma.^(3/2)) ./ ...
    (1 - gamma.^2).^(1/2) + 0.83 - 1.76 .* gamma) .* ...
    (1 - (1 - gamma).*lambda);
