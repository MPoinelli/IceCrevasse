% Mattia Poinelli
% JPL, April 2017
% Visualization of Digitalised Features and Smoothing elaboration
% TO SAVE the smoothed version uncomment 'shapewrite' line
% REMEMBER TO COMMENT WHEN DONE

clear all, close all, clc

theta    = 90 - [0:1:180]; theta_rad = deg2rad(theta); % latitudes
phi      = [0:1:360];  phi_rad = deg2rad(phi);         % longitudes

% Load Raster & plot map
[raster] = imread('Europa_raster_coarse.jpg');

figure(1),hold on
set(figure(1),'units','pixels','position',[0,500,1000,500])
axesm ('eqdcylin', 'Frame', 'on', 'Grid', 'on','MapLatLim',[-75,75],'MapLonLim',[-360,0]);
framem on, gridm on, mlabel on, plabel on
[lon,lat] = meshgrid(phi,theta); 
geoshow(lat,lon,raster);

S = shaperead('Features.shp');% ,...

C = struct2cell(S);

[m,n] = size(C);

for i = 1 : n 
    X = C{3,i}; X = X(~isnan(X));
    Y = C{4,i}; Y = Y(~isnan(Y));
    [X, index] = unique(X); 
    xx = X(1): 0.1 : X(end);
    yy = spline(X, Y(index), xx);
    S(i).X = xx;
    S(i).Y = yy;
end

%% UNCOMMENT if the smoothed shapefile is desired
%shapewrite(S,'Features_smoothed.shp')
geoshow([S(:).Y],[S(:).X]+180,...
    'DisplayType','Point','Marker','.','Color','red')