function beta = beta_theta_theta (theta,m,h,l)
% Beta-function evaluation as function of:
% 
% theta: co-latitude
% m:     order (0 to 2)
% h & l: Love numbers

if abs(m) > 2 
    error ('Order too high! chose a value lower than 2!')
end

if m == 0
    
    beta = (3/4)*(3*h-10*l)*cos(2*theta) + (3/4)*(h-2*l);
    
elseif m == 1
    
    beta = (3/2)*(3*h-10*l)*sin(2*theta);

elseif m == 2

    beta = -(3/2)*(3*h-10*l)*cos(2*theta) + (9/2)*(h-2*l);
    
end

end
