% Ir. Mattia Poinelli
% Verona, October 2017
%
% Calculation of shape_file length
% 
% KEY OUTPUT:
%
% vector LENGTH containing the lengths of the crevasses

%% Load parameters & shape file
run ../_physical_parameters/Europa_physics.m
run ../_physical_parameters/Europa_LEFM_parameters.m

S = shaperead('../_target_features/Features_smoothed.shp');

%% Lenght loop
for q = 1:26
    
    Length_crevasse = 0;
    % avoid NaN values, xx = longitudes and yy = latitudes
    
    xx = S(q).X; xx = xx(~isnan(xx));
    yy = S(q).Y; yy = yy(~isnan(yy));
    END = length(xx);
    
    for i = 1:END
        
        %% Rotation matrices
        if i == length(xx)
            % last node of the grid has the same charateristic of the
            % previous node to avoid Matlab unconsistencies
            ROT(:,:,i) =  ROT(:,:,i-1);
            
            Tang_ve(1,i) = Tang_ve(1,i-1);
            Tang_ve(2,i) = Tang_ve(2,i-1);
            Norm_ve(1,i) = Norm_ve(1,i-1);
            Norm_ve(2,i) = Norm_ve(2,i-1);
            
            angular_distance (i) = angular_distance(i-1);

        elseif yy(i+1) > yy(i) % 'ascending' feature
            % evaluation of useful angles a,b
            b_ang(i) = deg2rad(xx(i+1) - xx(i));
            a_ang(i) = deg2rad(yy(i+1) - yy(i));
            
            % Side_Angle_Side problem, A is the rotation angle, and the
            % last output c gives the length of the segment
            [A,~,c] = Side_Angle_Side (a_ang(i),b_ang(i));
            
            ROT(:,:,i) = [cos(A),-sin(A);
                sin(A),cos(A)]; % counter clock-wise rotation (positive angle)
            
            % Tangent vector
            Tang_ve(1,i) = ROT(1,1,i);
            Tang_ve(2,i) = ROT(2,1,i);
            
            % Normal vector
            Norm_ve(1,i) = ROT(1,2,i);
            Norm_ve(2,i) = ROT(2,2,i);
            
            angular_distance (i) = c;
            
        elseif yy(i+1) < yy(i)% 'descending' feature
            % evaluation of useful angles a,b
            b_ang(i) = deg2rad(xx(i+1) - xx(i));
            a_ang(i) = abs(deg2rad(yy(i+1) - yy(i)));
            
            % Side_Angle_Side problem, A is the rotation angle, and the
            % last output c gives the length of the segment
            [A,~,c] = Side_Angle_Side (a_ang(i),b_ang(i));
            
            ROT(:,:,i) = [cos(-A),-sin(-A);
                sin(-A),cos(-A)]; % clock-wise rotation (negative angle)
            
            % Tangent vector
            Tang_ve(1,i) = ROT(1,1,i);
            Tang_ve(2,i) = ROT(2,1,i);
            
            % Normal vector
            Norm_ve(1,i) = ROT(1,2,i);
            Norm_ve(2,i) = ROT(2,2,i);
            
            angular_distance (i) = c;
        end
        
        Length_crevasse = Length_crevasse + angular_distance (i)*R;
        
    end
    Mean_segment (q)= mean(angular_distance*R)/1000;
    LENGTH (q) = Length_crevasse/1000;
    clear Length_crevasse xx yy
     
    Name = strcat(num2str(S(q).Name));
    DONE = strcat(num2str(S(q).Name),'''s calculated length is :',num2str(ceil(LENGTH(q))),' km');
    disp(DONE)
end