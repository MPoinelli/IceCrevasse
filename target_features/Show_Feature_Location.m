% Mattia Poinelli
% JPL, April 2019
% Visualization of Digitalised Features

%clear all, close all, clc

theta    = 90 - [0:1:180]; theta_rad = deg2rad(theta); % latitudes
phi      = [0:1:360];  phi_rad = deg2rad(phi);         % longitudes

% Load Raster & plot map
[raster] = imread('Europa_raster_coarse.jpg');

figure(1),hold on
set(figure(1),'units','pixels','position',[0,500,1000,500])
%axesm ('eqdcylin', 'Frame', 'on', 'Grid', 'on','MapLatLim',[-75,75],'MapLonLim',[-360,0]);
%framem on, gridm on, mlabel on, plabel on
[lon,lat] = meshgrid(phi,theta); 
geoshow(lat,lon+180,raster);

S = shaperead('Cycloids_lon_lat.shp');% ,...
C = struct2cell(S);
% 
 [m,n] = size(C);
% % Number of Nodes
% % 100, 1000 or 10000
  Nodes = 100;
% 
% % Generation of multiple nodes
for i = 1 : n
    X = C{3,i}; X = X(~isnan(X));
    Y = C{4,i}; Y = Y(~isnan(Y));
    [X, index] = unique(X);
    
    
    [pt] = interparc(Nodes,X,Y,'linear');
    
    %     xx = X(1): 0.01 : X(end);
    %     yy = spline(X, Y(index), xx);
    S(i).X = pt(:,1)';
    S(i).Y = pt(:,2)';
    
    % Cycloid numbering problem
    if i == 2
        xx_save = pt(:,1)';
        yy_save = pt(:,2)';
    elseif i == 5
        S(5).X = xx_save;
        S(5).Y = yy_save;
        S(2).X = pt(:,1)';
        S(2).Y = pt(:,2)';
    elseif i == 6
        S(6).X = xx_save;
        S(6).Y = yy_save;
        S(5).X = pt(:,1)';
        S(5).Y = pt(:,2)';
    end
    clear pt
end

geoshow([S(:).Y],[S(:).X],...
    'DisplayType','Point','Marker','o','Color','red')