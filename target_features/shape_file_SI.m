% Mattia Poinelli
% JPL, May 2019
% Simple routine to produce the shape files included in the SI

clear all, close all, clc

S = shaperead('Cycloids_smoothed10000.shp');% ,...

S_new(1) = S(1);
S_new(2) = S(2);
S_new(3) = S(4);
S_new(4) = S(6);

theta    = 90 - [0:1:180]; theta_rad = deg2rad(theta); % latitudes
phi      = [0:1:360];  phi_rad = deg2rad(phi);         % longitudes

% Load Raster & plot map
[raster] = imread('Europa_raster_coarse.jpg');

figure(1),hold on
set(figure(1),'units','pixels','position',[0,500,1000,500])
%axesm ('eqdcylin', 'Frame', 'on', 'GrCRAid', 'on','MapLatLim',[-75,75],'MapLonLim',[-360,0]);
%framem on, gridm on, mlabel on, plabel on
axesm eqdcylin
[lon,lat] = meshgrid(phi,theta);
geoshow(lat,lon,raster);
geoshow([S_new(:).Y],[S_new(:).X],...
    'DisplayType','Point','Marker','o','Color','red')

% Positive West 0-360
S_new(1).X(:) = S_new(1).X(:)+180;
S_new(2).X(:) = S_new(2).X(:)+180;
S_new(3).X(:) = S_new(3).X(:)+180;
S_new(4).X(:) = S_new(4).X(:)+180;

for i= 1:4
    textm ([S_new(i).Y(1)],[S_new(i).X(1)],num2str(i))
end