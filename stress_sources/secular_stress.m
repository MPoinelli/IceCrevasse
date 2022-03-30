function [sigma_theta,sigma_phi,tau] = secular_stress (theta,phi,n,t,T_ns,g,R,mu,eta,h,l)
% Analytical Secular Stress
% Poinelli Mattia, JPL, March 2017
%
% Analytical formulation of diurnal stress that includes effects of
% non-synchronous rotation
% Equations from Jara-Orue, 2011
%
% INPUT
%
% theta, phi: COlatitude & longitude [rad]
% n:          mean motion            [rad/sec]
% t:          time after pericenter passage [sec]
% T_ns:       period of NS rotation  [sec]
% g:          gravity at the surface [m/s2]
% R:          radius                 [m]
% mu:         ridigidy               [Pa]
% eta:        viscosity              [Pa s]
% h,l:        Love parameters
% 
% OUTPUT
% 
% sigma_theta: colatitudinal stress [Pa]
% sigma_phi:   longitudinal stress  [Pa]
% tau:         shear stress         [Pa]

% Declare global variables
global Delta Z

% Define global variables
Delta = T_ns/(4*pi*(eta/mu));
Z = 0.5 * (n^2*R*mu)/(g*sqrt(1+Delta^2));

% Colatitudinal Direction
sigma_theta = Z * alpha_theta_theta (theta,2,h,l) * ...
    cos(2*phi + 4 * pi * t /T_ns + atan(Delta));

% Longitudinal Direction
sigma_phi = Z * alpha_phi_phi (theta,2,h,l) * ...
     cos(2*phi + 4 * pi * t /T_ns + atan(Delta));
 
% Tangential Direction
tau = - Z * alpha_theta_phi (theta,2,h,l) * ...
    sin(2*phi + 4 * pi * t /T_ns + atan(Delta));

end
