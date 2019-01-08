function beta = beta_theta_phi (theta,m,h,l)
% Beta-function evaluation as function of:
% 
% theta: co-latitude
% m:     order (0 to 2)
% h & l: Love numbers

if abs(m) > 2 
    error ('Order too high! chose a value lower than 2!')
elseif abs(m) == 0
    error ('Order zero not required! chose a value lower than 2!')
end
    
if m == 1
    
    beta = 3*l*sin(theta);

elseif m == 2

    beta = 3*l*cos(theta);
    
end

end