% Mattia Poinelli
% JPL, May 2017

% VISUALISATION OF GLOBAL MAP for different purposes
% Currently loading: mean global depths.

clear all,  clc, close all

run Europa_physics.m

% Load Raster & plot map
figure('Units','inches','Position',[0 0 7 4.5],'PaperPositionMode','auto')

%depth (1,:) =  max(max(depth)); depth (end,: ) =  max(max(depth));
load('Critical_depths_file.mat')
depth = depth(:,:,1:end-2);
Mean_depth = mean(depth,3);

% generation of grid
phi      = [0:20:360]; phi_rad = deg2rad(phi);          % longitudes
theta    = 90 - [0:20:180]; theta_rad = deg2rad(theta); % latitudes

% c = colorbar;
% c.Label.String = 'Average Critical Depth [m]';c.FontSize = 15;
% c.Limits = [0 70];
% c.FontName = 'latex';

xt = get (gca,'XTick');
yt = get (gca,'YTick');
 
for k=1:numel(xt);
xt1{k}=sprintf('%d°',xt(k));
end
 
for k=1:numel(yt);
yt1{k}=sprintf('%d°',yt(k));
end

% Load Raster & plot map
[raster] = imread('Europa_raster_coarse.jpg');
S = shaperead('Features_smoothed.shp');

%depth (1,:) =  max(max(depth)); depth (end,: ) =  max(max(depth));

%Tit = title ('Average Critical Depth');
% set(Tit,'Fontsize',8)
% set(Tit,'Interpreter','latex'))
axesm('MapProjection','eqdcylin','FontName','times','FontSize',11,...
    'MapLatLimit',[-90 90],'MapLonLimit', [0 360],...
    'LabelFormat','none',...
    'MLabelLocation',30,'PLabelLocation',30)

framem on, plabel on, mlabel (-85)
[lon,lat] = meshgrid(phi,theta);

contourfm (theta,phi,Mean_depth,'Linestyle','none');
Max_mean_height = 871.6667;

cs = colorbar;
cs.FontName = 'Times';
cs.Location = 'Northoutside'
%c.Position = [0.9165    0.1577    0.0310    0.7100]
%c.Label.String = 'Average Critical Depth [m]';c.FontSize = 12;
cs.Limits = [0 max(max(Mean_depth))];
cs.Label.String = 'Depth of Surface Crevasses';
cs.TickLabels = {'0 m','10 m','20 m','30 m','40 m','50 m','60 m','70 m'};

cb = colorbar;
cb.FontName = 'Times';
cb.Location = 'SouthOutside'
%c.Position = [0.9165    0.1577    0.0310    0.7100]
%c.Label.String = 'Average Critical Depth [m]';c.FontSize = 12;
%cb.Limits = [0 87.16667];
cb.Limits = [0 max(max(Mean_depth))];
cb.Ticks = linspace (0, max(max(Mean_depth)),10); cb.Ticks = cb.Ticks(1:9);
 cb.TickLabels = {'0 m','100 m','200 m','300 m','400 m','500 m','600 m','700 m','800 m'};
cb.Label.String = 'Height of Bottom Crevasses';

tightmap, gridm on



geoshow([S([1,2,5,6,12,15,16,21,24]).Y],[S([1,2,5,6,12,15,16,21,24]).X]+180,...
    'Color','white','Linewidth',2)

NUMBERS = [1 2 5 6 4 9 8 7 3];
j=0;
for i =[1,2,5,6,12,15,16,21,24]
    
    j = j+1;
    L = length(S(i).X(:));
    XX(j) = S(i).X(ceil(L/2));
    YY(j) = S(i).Y(ceil(L/2));
    textm(YY(j)+11,XX(j)+180,num2str(NUMBERS(j)),'color','white','FontSize',17)

end
tightmap
% textm(YY(1),XX(1)+180,'1','white')
% textm(YY(2),XX(2)+180,'2','white')
% textm(YY(3),XX(3)+180,'3','white')
% textm(YY(4),XX(4)+180,'4','white')
% textm(YY(5),XX(5)+180,'5','white')
% textm(YY(6),XX(6)+180,'6','white')
% textm(YY(7),XX(7)+180,'7','white')
% textm(YY(8),XX(8)+180,'8','white')
% textm(YY(9),XX(9)+180,'9','color','white')
