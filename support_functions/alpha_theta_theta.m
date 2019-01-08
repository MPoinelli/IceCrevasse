function alpha = alpha_theta_theta (theta,m,h,l)
% Alpha-function evaluation as function of:
% 
% theta: co-latitude
% m:     order (0 to 2)
% h & l: Love numbers

if m ~= 2 
    error ('Only order 2 is acceptable!')
    
end

alpha = -(3/2)* (3*h-10*l)* cos(2*theta) + (9/2)*(h-2*l);