function [sigma_theta,sigma_phi,tau] = diurnal_stress (theta,phi,n,e,w,t,epsilon,g,R,mu,eta,h,l)
% Analytical Diurnal Stress
% Poinelli Mattia, JPL, March 2017
%
% Analystical formulation of diurnal stress that includes effects of
% non-zero eccentricity and non-zero obliquity
% Equations  from Jara-Orue, 2011
%
% INPUT
%
% theta, phi: COlatitude & longitude [rad]
% n:          mean motion            [rad/sec]
% e:          eccentricity           [~]
% w:          Argument of pericenter [rad]
% t:          time after pericenter passage [sec]
% epsilon:    obliquity              [rad]
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

% Constant
Lambda = (mu/eta)/n;
Z = 0.5 * (n^2*R*mu)/(g*sqrt(1+Lambda^2));

% Colatitudinal Direction
sigma_theta = Z * (-6*e*beta_theta_theta(theta,0,h,l)*cos(n*t+atan(Lambda)) + ...
    e*beta_theta_theta(theta,0,h,l)*(4*sin(2*phi)* sin (n*t + atan(Lambda)) + ...
    3*cos(2*phi)*cos(n*t+atan(Lambda))) + ...
    4*cos(epsilon)*sin(epsilon)* beta_theta_theta(theta,1,h,l) * ...
    (cos(phi)*sin(w + n*t + atan(Lambda))));

% Longitudinal Direction
sigma_phi = Z * (-6*e*beta_phi_phi(theta,0,h,l)*cos(n*t+atan(Lambda)) + ...
    e*beta_phi_phi(theta,0,h,l)*(4*sin(2*phi)* sin (n*t + atan(Lambda)) + ...
    3*cos(2*phi)*cos(n*t+atan(Lambda))) + ...
    4*cos(epsilon)*sin(epsilon)* beta_phi_phi(theta,1,h,l) * ...
    (cos(phi)*sin(w + n*t + atan(Lambda))));

% Tangential Direction
tau = Z* (2*e*beta_theta_phi(theta,2,h,l) * ...
    (4* cos(2*phi)* sin(n*t + atan(Lambda)) - ...
    3* sin(2*phi)*cos(n*t+atan(Lambda))) + ...
    4 * cos(epsilon)*sin(epsilon)* beta_theta_phi(theta,1,h,l)* ...
    sin(phi)*sin(w + n*t + atan(Lambda)));

end