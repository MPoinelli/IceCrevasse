% Physical parameters of Europa
% from Jara-Orue&Vermeersen, 2011 

G = 6.67408e-11; % [m3 kg-1 s-2] gravitational constant

%% Orbital parameters
n = deg2rad(101.37472)/(60*60*24); % [rad/s] mean motion
T_europa = 2*pi/n;                 % [sec]   orbital period
e = 0.0094;                        % [~]     eccentricity
w = deg2rad(345);                  % [rad]   argument of pericenter 

%% Physical parameters
T_ns_years = 10e7;             % [years] period of Non-Synchronous rotation
T_ns = T_ns_years*365*24*3600; % [sec]   same as above
R = 1562000;                   % [m]     radius of Europa
epsilon = deg2rad(0.1);        % [rad]   obliquity [Bills, 2005]
g = 1.315;                     % [m/s2]  gravity at the surface
Mass = 4.7998e22;              % [Kg]    mass of Europa

%% Rheology of the Crust
Bulk_Modulus = 9.3e9;    % [Pa]  Bulk modulus
mu = 3.487e9;            % [Pa]  rigidity 
E = 9*Bulk_Modulus*mu/...
    (3*Bulk_Modulus+mu); % [Pa]  Youngs modulus
eta = 10^21;             % [Pas] Viscosity
nu = E/(2*mu)-1;         % [~]   Poisson ratio

h_d = 1.151;           % Love number diurnal
l_d = 3.07996*10^(-1); % Love number diurnal
h_s = 1.85155;         % Love number secular
l_s = 4.95366*10^(-1); % Love number secular
