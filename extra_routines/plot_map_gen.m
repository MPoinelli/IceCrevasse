function plot_map(lambda, theta, grid, cmin,cmax)
% Plot gridded data on Earth's surface.
%
%.i/ lambda: p x 1 vector containing the longitudes (in radian).
%.i/ theta: q x 1 vector containing the co-latitudes (in radian).
%.i/ grid: q x p matrix containing the gridded values.
%.i-optional/ cmin,cmax: sets the color limits to specified minimum and
% maximum values respectively


% check dimensions
cols = numel(lambda);
if ~isequal(size(lambda),[1,cols])
   lambda = reshape(lambda,1,cols);
end

rows = numel(theta);
if ~isequal(size(theta),[1,rows])
   theta = reshape(theta,1,rows);
end

if size(grid,1) ~= rows || size(grid,2) ~= cols
    error('check dimension!')
end

[lon,lat] = meshgrid(lambda*180/pi,90 - theta*180/pi);
h = figure;set(gca,'FontSize',14,'FontName','Arial');axesm eckert4;framem; gridm;axis off    
%hcb = colorbar('horiz');hold on;
%load('coastlines');
pcolorm(lat,lon,grid);
%plotm(coastlat, coastlon,'k')


if nargin > 3
set(get(h,'currentaxes'),'clim',[cmin,cmax])
end
