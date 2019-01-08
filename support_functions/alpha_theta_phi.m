function alpha = alpha_theta_phi (theta,m,h,l)
% Alpha-function evaluation as function of:
% 
% theta: co-latitude
% m:     order (0 to 2)
% h & l: Love numbers

if m ~= 2 
    error ('Only order 2 is acceptable!')
    
end

alpha = 3*l*cos(theta);

